"""Scenario analysis to inform cross-sector benchmark update, July 2023."""
import os

import numpy
import pandas

# directory containing scenario data
_PROJ_DIR = "C:/Users/Ginger.Kowal.CORPCDP/Documents/Scenario review"

# output directory
_OUT_DIR = os.path.join(
    os.path.dirname(__file__), 'output')

# GWP100 conversion values from AR5 report
_N2O_GWP100_AR5 = 265
_CH4_GWP100_AR5 = 28

# GW100 conversion values from AR6
_CH4FOSS_GWP100_AR6 = 29.8  # GWP100 for fossil methane
_N2O_GWP100_AR6 = 273
_HFC_GWP100_AR6 = 1530
_PFC_GWP100_AR6 = 7380
_SF6_GWP100_AR6 = 25200

# cumulative emissions from FLAG, 2020-2050 (GtCO2e)
_FLAG_2020_50_CO2E = -99.54

# conversion factor from kt to Mt
_KT_to_MT = 0.001

# conversion factor from Gt to Mt
_GT_to_MT = 1000

# medium concern, yearly energy from bioenergy in any year between 2010-2050
# (EJ/year)
_MED_BIO = 100.

# maximum yearly sequestration via af-/reforestation in any year between
# 2010-2050 (Gt CO2/year)
_MAX_AFOLU = 3.6

# maximum cumulative CCS between 2010 and 2050 (Gt)
_MAX_CCS = 214

# gross energy and industrial process CO2 emissions in 2020 (Mt CO2)
# (Forster et al 2023)
_2020_CO2 = 35289.751


def read_ar6_data():
    """Read AR6 scenario data from file.

    Args:
        None

    Returns:
        a tuple, (ar6_key, ar6_scen) containing:
        ar6_key (Pandas dataframe): dataframe containing scenario metadata
        ar6_scen (Pandas dataframe): dataframe containing scenario data

    """
    key_path = os.path.join(
        _PROJ_DIR, 'IPCC_AR6',
        'AR6_Scenarios_Database_metadata_indicators_v1.1.xlsx')
    ar6_key = pandas.read_excel(
        key_path, sheet_name='meta_Ch3vetted_withclimate')
    ar6_key['scen_id'] = ar6_key['Model'] + ' ' + ar6_key['Scenario']

    scen_path = os.path.join(
        _PROJ_DIR, 'IPCC_AR6/AR6_Scenarios_Database_World_v1.1.csv')
    ar6_scen = pandas.read_csv(scen_path)
    ar6_scen['scen_id'] = ar6_scen['Model'] + ' ' + ar6_scen['Scenario']

    return ar6_key, ar6_scen


def sustainability_filters(scen_id_list, emissions_df):
    """Filter scenarios according to sustainability thresholds.

    Filter a set of scenarios according to 2050 deployment of biofuels, CCS,
    and af-/reforestation. Deployment in 2050 is calculated as the average of
    deployment in 2040 and 2060.

    Args:
        scen_id_list (list of strings): list of scenario ids that should be
            filtered. For example, this could be the full list of scenarios
            in the AR6 database, or the subset of C1 scenarios
        emissions_df (Pandas dataframe): dataframe containing data for
            emissions and sequestration, used to identify scenarios meeting
            filtering criteria

    Returns:
        a list of strings that is a subset of `scen_id_list`, giving the
            scenarios that meet the given filter criteria

    """
    test_col = [str(idx) for idx in list(range(2010, 2051))]
    biom_df = emissions_df.loc[
        emissions_df['Variable'] == 'Primary Energy|Biomass']
    rem_biom = set(
        biom_df.loc[(biom_df[test_col] > _MED_BIO).any(axis=1)]['scen_id'])

    lu_var = 'Carbon Sequestration|Land Use|Afforestation'
    lu_df = emissions_df.loc[emissions_df['Variable'] == lu_var]
    rem_lu = set(
        lu_df.loc[
            (lu_df[test_col] > (_MAX_AFOLU * _GT_to_MT)).any(axis=1)][
            'scen_id'])

    ccs_df = emissions_df.loc[
        emissions_df['Variable'] == 'Carbon Sequestration|CCS']
    rem_ccs = set(
        ccs_df.loc[ccs_df[test_col].interpolate(
            axis=1).sum(axis=1) > (_MAX_CCS * _GT_to_MT)]['scen_id'])

    rem_afolu = flag_filter(emissions_df)

    # filter scenarios with >2.3 Mt CDR in 2020
    ccs_var_list = [
        'Carbon Sequestration|CCS|Biomass',
        'Carbon Sequestration|Direct Air Capture',
        'Carbon Sequestration|Enhanced Weathering']
    cdr_df = emissions_df.loc[emissions_df['Variable'].isin(ccs_var_list)]
    cdr_sum_df = cdr_df.groupby('scen_id').sum()
    cdr_sum_df.reset_index(inplace=True)
    rem_cdr = set(cdr_sum_df.loc[cdr_sum_df['2020'] > _2020_CDR]['scen_id'])

    rem_ids = rem_biom.union(rem_lu).union(rem_ccs).union(rem_afolu).union(
            rem_cdr)

    filtered_ids = set(scen_id_list).difference(rem_ids)
    return filtered_ids


def fill_EIP_emissions(em_df, id_list):
    """Calculate net energy & industrial process emissions.

    Not all scenarios contain the variable "Emissions|CO2|Energy and Industrial
    Processes". For those that don't, calculate it as the sum of energy, and
    industrial process CO2 emissions.

    Args:
        em_df (Pandas dataframe): dataframe containing emissions data
        id_list (list): complete list of scenario ids to fill

    Returns:
        a Pandas dataframe containing all the values in `em_df`, plus
            calculated values for Energy and Industrial Process emissions
            for those scenarios that didn't originally report it

    """
    year_col = [col for col in em_df if col.startswith('2')]
    summary_cols = ['scen_id', 'Variable'] + year_col
    co2_var = 'Emissions|CO2|Energy and Industrial Processes'
    co2_em = em_df.loc[em_df['Variable'] == co2_var][summary_cols]
    em_df.reset_index(inplace=True)

    # not all models report values for this variable: identify
    # those that don't
    estimate_mod = set(id_list).difference(
        set(em_df.loc[em_df['Variable'] == co2_var]['scen_id']))
    print("**note: estimating EIP emissions for {} scenarios".format(
        len(estimate_mod)))

    # for those, calculate the variable as the sum of emissions from Energy,
    # and from Industrial processes
    r_idx = len(em_df)
    for sid in estimate_mod:
        sum_rows = em_df.loc[
            (em_df['scen_id'] == sid) &
            (em_df['Variable'].isin(
                ['Emissions|CO2|Energy',
                'Emissions|CO2|Industrial Processes']))]
        sum_vals = sum_rows[year_col].sum()
        em_df.loc[r_idx] = sum_vals
        em_df.loc[r_idx, 'scen_id'] = sid
        em_df.loc[r_idx, 'Variable'] = co2_var
        r_idx = r_idx + 1

    return em_df


def calc_gross_eip_co2(em_df, year_col):
    """Calculate gross CO2 emissions from EIP.

    Gross CO2 emissions from energy and industrial prcoesses are calculated
    as `Emissions|CO2|Energy and Industrial Processes` +
    `Carbon Sequestration|CCS|Biomass` +
    `Carbon Sequestration|Direct Air Capture` +
    `Carbon Sequestration|Enhanced Weathering`

    Args:
        em_df (pandas dataframe): dataframe containing emissions variables
        year_col (list): list of columns giving yearly values

    Returns:
        dataframe containing gross EIP CO2 emissions

    """
    ccs_var_list = [
        'Carbon Sequestration|CCS|Biomass',
        'Carbon Sequestration|Direct Air Capture',
        'Carbon Sequestration|Enhanced Weathering']
    # check for negative values in CCS. if any values are negative,
    # change the sign
    ccs_df = em_df.loc[em_df['Variable'].isin(ccs_var_list)]
    if (ccs_df[year_col] < 0).any().any():
        for year in year_col:
            ccs_df.loc[ccs_df[year] < 0, year] = -(
                ccs_df.loc[ccs_df[year] < 0, year])
    eip_var = 'Emissions|CO2|Energy and Industrial Processes'
    eip_df = em_df.loc[em_df['Variable'] == eip_var]
    sum_df = pandas.concat([ccs_df, eip_df])
    sum_cols = year_col + ['scen_id']

    gross_eip_co2_df = sum_df[sum_cols].groupby('scen_id').sum()
    gross_eip_co2_df.reset_index(inplace=True)
    gross_eip_co2_df['Variable'] = 'Emissions|CO2|Energy and Industrial Processes|Gross'
    return gross_eip_co2_df


def calc_eip_n2o(em_df, year_col):
    """Calculate N2O emissions from energy and industrial processes.

    Energy-related (i.e., non-AFOLU) N2O emissions are calculated as:
    'Emissions|N2O|Energy' + 'Emissions|N2O|Industrial Processes' +
    'Emissions|N2O|Other' + 'Emissions|N2O|Waste'.

    Args:
        em_df (pandas dataframe): dataframe containing emissions variables
        year_col (list): list of columns giving yearly values

    Returns:
        A pandas dataframe containing energy-related N2O emissions for each
            scenario, in kt N2O per year

    """
    n2o_var_list = [
        'Emissions|N2O|Energy', 'Emissions|N2O|Industrial Processes',
        'Emissions|N2O|Other', 'Emissions|N2O|Waste']
    sum_cols = year_col + ['scen_id']
    n2o_df = em_df.loc[em_df['Variable'].isin(n2o_var_list)]
    eip_n2o_df = n2o_df[sum_cols].groupby('scen_id').sum()
    return eip_n2o_df


def calc_eip_ch4(em_df, year_col):
    """Calculate CH4 emissions from energy and industrial processes.

    Energy-related (i.e., non-AFOLU) CH4 emissions are calculated as:
    'Emissions|CH4|Energy' + 'Emissions|CH4|Industrial Processes' +
    'Emissions|CH4|Other' + 'Emissions|N2O|Waste'.

    Args:
        em_df (pandas dataframe): dataframe containing emissions variables
        year_col (list): list of columns giving yearly values

    Returns:
        A pandas dataframe containing energy and industrial process CH4
            emissions for each scenario, in Mt CH4 per year

    """
    ch4_var_list = [
        'Emissions|CH4|Energy', 'Emissions|CH4|Industrial Processes',
        'Emissions|CH4|Other', 'Emissions|CH4|Waste']
    sum_cols = year_col + ['scen_id']
    ch4_df = em_df.loc[em_df['Variable'].isin(ch4_var_list)]
    eip_ch4_df = ch4_df[sum_cols].groupby('scen_id').sum()
    return eip_ch4_df


def flag_filter(emissions_df):
    """Filter scenarios for compatibility with SBTi FLAG pathway.

    The SBTi FLAG pathway uses mitigation potentials from Roe et al (2019) and
    a baseline emissions value (from Roe et al) to calculate a 1.5C-compatible
    emissions pathway for FLAG. Use this pathway to identify scenarios where
    emissions from the land sector are smaller than this in terms of cumulative
    CO2e between 2020 and 2050.  CH4 and N2O are converted to CO2eq using
    GWP-100 values from AR5 (comparable to the FLAG pathway).

    Args:
        emissions_df (Pandas dataframe): dataframe containing data for
            emissions and sequestration, used to identify scenarios meeting
            filtering criteria

    Returns:
        a list of strings, giving the scenarios that should be removed
            according to compatibility with the FLAG pathway

    """
    test_col = [str(idx) for idx in list(range(2020, 2051))]
    # calculate total AFOLU CO2e
    sum_cols = test_col + ['scen_id']
    co2_df = emissions_df.loc[
        emissions_df['Variable'] == 'Emissions|CO2|AFOLU']
    ch4_df = emissions_df.loc[
        emissions_df['Variable'] == 'Emissions|CH4|AFOLU']
    ch4_co2eq = ch4_df[sum_cols] * _CH4_GWP100_AR5
    n2o_df = emissions_df.loc[
        emissions_df['Variable'] == 'Emissions|N2O|AFOLU'][
            sum_cols].groupby('scen_id').sum()
    n2o_co2eq = n2o_df * _N2O_GWP100_AR5 * _KT_to_MT
    n2o_co2eq.reset_index(inplace=True)

    co2e_df = pandas.concat([co2_df, ch4_df, n2o_df]).groupby('scen_id').sum()
    co2e_df.replace(0, numpy.nan, inplace=True)
    co2e_df.reset_index(inplace=True)

    afolu_em = co2e_df[test_col].interpolate(axis=1).sum(axis=1)
    afolu_limit = _FLAG_2020_50_CO2E  * _GT_to_MT
    rem_afolu = set(
        co2e_df.loc[co2e_df[test_col].interpolate(
            axis=1).sum(axis=1) < afolu_limit]['scen_id'])

    return rem_afolu


def summarize_c1_key_var():
    """Summarize key variables in AR6 C1 scenarios."""
    ar6_key, ar6_scen = read_ar6_data()
    year_col = [col for col in ar6_scen if col.startswith('2')]

    c1_scen = ar6_key.loc[ar6_key['Category'] == 'C1']['scen_id']
    c1_em = ar6_scen.loc[ar6_scen['scen_id'].isin(c1_scen)]
    c1_em.set_index('scen_id', inplace=True)
    summary_cols = ['2050', 'Variable']

    # Final energy demand in 2050
    final_energy = c1_em.loc[c1_em['Variable'] == 'Final Energy']['2050']

    # Maximum yearly primary energy from bioenergy, 2010-2050
    test_col = [str(idx) for idx in list(range(2010, 2051))]
    biom_df = c1_em.loc[c1_em['Variable'] == 'Primary Energy|Biomass']
    max_biom = biom_df[test_col].max(axis=1)

    # Cumulative CCS, 2010-2050
    ccs_df = c1_em.loc[c1_em['Variable'] == 'Carbon Sequestration|CCS']
    cum_ccs = ccs_df[test_col].interpolate(axis=1).sum(axis=1)

    # Total atmospheric CDR in 2050
    cdr_var_list = [
        'Carbon Sequestration|CCS|Biomass',
        'Carbon Sequestration|Direct Air Capture',
        'Carbon Sequestration|Enhanced Weathering']
    cdr_df = c1_em.loc[
        c1_em['Variable'].isin(cdr_var_list)][summary_cols]
    cdr_sum = cdr_df.groupby('scen_id').sum()['2050']

    # Share of primary energy from renewables in 2050
    c1_prien_2050 = c1_em.loc[c1_em['Variable'] == 'Primary Energy']['2050']
    c1_renen_2050 = c1_em.loc[
        c1_em['Variable'] == 'Primary Energy|Renewables (incl. Biomass)'][
            '2050']
    c1_en_df = pandas.DataFrame({
        'Primary Energy': c1_prien_2050,
        'Primary Energy|Renewables': c1_renen_2050})
    ren_share_2050 = (c1_en_df['Primary Energy|Renewables'] /
        c1_en_df['Primary Energy'])

    # Cumulative gross fossil CO2, 2020-2050
    em_col = [str(idx) for idx in list(range(2020, 2051))]
    ar6_filled_em = fill_EIP_emissions(ar6_scen, c1_scen)
    ar6_gross_em = calc_gross_eip_co2(ar6_filled_em, year_col)
    c1_gross_em = ar6_gross_em.loc[ar6_gross_em['scen_id'].isin(c1_scen)]
    c1_gross_em.set_index('scen_id', inplace=True)
    c1_gross_em.replace(0, numpy.nan, inplace=True)
    cum_gross_co2 = c1_gross_em[em_col].interpolate(axis=1).sum(axis=1)

    key_var_vals = pandas.DataFrame({
        'Primary energy 2050': final_energy,
        'Max yearly bioenergy': max_biom,
        'Cumulative CCS 2010-2050': cum_ccs,
        'CDR 2050': cdr_sum,
        'Ren share 2050': ren_share_2050,
        'Cumulative gross fossil CO2 2020-2050': cum_gross_co2,
        })
    key_var_25p = key_var_vals.quantile(q=0.25)
    key_var_75p = key_var_vals.quantile(q=0.75)
    key_var_df = pandas.DataFrame({
        'C1 25 perc': key_var_25p,
        'C1 75 perc': key_var_75p,
        })
    key_var_df.to_csv(
        os.path.join(_OUT_DIR, 'AR6_C1_key_var_summary.csv'))


def calc_afolu_co2(emissions_df):
    id_df = emissions_df.set_index('scen_id')
    test_col = [str(idx) for idx in list(range(2020, 2051))]
    # calculate total AFOLU CO2e
    sum_cols = test_col + ['scen_id']
    co2_df = id_df.loc[id_df['Variable'] == 'Emissions|CO2|AFOLU']
    co2_df.replace(0, numpy.nan, inplace=True)
    afolu_co2 = co2_df[test_col].interpolate(axis=1).sum(axis=1)
    afolu_co2_df = pandas.DataFrame(afolu_co2)
    afolu_co2_df.reset_index(inplace=True)
    return afolu_co2_df

    
def summarize_filtered_key_var():
    """Summarize key variables in the final scenario set."""
    ar6_key, ar6_scen = read_ar6_data()
    year_col = [col for col in ar6_scen if col.startswith('2')]

    c1_scen = ar6_key.loc[ar6_key['Category'] == 'C1']['scen_id']
    c1_filtered = sustainability_filters(c1_scen, ar6_scen)
    filtered_em = ar6_scen.loc[ar6_scen['scen_id'].isin(c1_filtered)]
    filtered_em.set_index('scen_id', inplace=True)
    summary_cols = ['2050', 'Variable']

    # Net CO2e AFOLU emissions in 2050
    afolu_co2e = calc_afolu_co2e(filtered_em, _CH4_GWP100_AR5, _N2O_GWP100_AR5)

    # Maximum yearly primary energy from bioenergy, 2010-2050
    test_col = [str(idx) for idx in list(range(2010, 2051))]
    biom_df = filtered_em.loc[
        filtered_em['Variable'] == 'Primary Energy|Biomass']
    max_biom = biom_df[test_col].max(axis=1)

    # final energy demand in 2050
    fin_energy_2050 = filtered_em.loc[
        filtered_em['Variable'] == 'Final Energy']['2050']

    # Cumulative CCS, 2010-2050
    ccs_df = filtered_em.loc[
        filtered_em['Variable'] == 'Carbon Sequestration|CCS']
    cum_ccs = ccs_df[test_col].interpolate(axis=1).sum(axis=1)

    # Cumulative gross fossil CO2, 2020-2050
    em_col = [str(idx) for idx in list(range(2020, 2051))]
    ar6_filled_em = fill_EIP_emissions(ar6_scen, c1_scen)
    ar6_gross_em = calc_gross_eip_co2(ar6_filled_em, year_col)
    filt_gross_em = ar6_gross_em.loc[ar6_gross_em['scen_id'].isin(c1_filtered)]
    filt_gross_em.set_index('scen_id', inplace=True)
    filt_gross_em.replace(0, numpy.nan, inplace=True)
    cum_gross_co2 = filt_gross_em[em_col].interpolate(axis=1).sum(axis=1)

    # Total atmospheric CDR in 2050
    cdr_var_list = [
        'Carbon Sequestration|CCS|Biomass',
        'Carbon Sequestration|Direct Air Capture',
        'Carbon Sequestration|Enhanced Weathering']
    cdr_df = filtered_em.loc[
        filtered_em['Variable'].isin(cdr_var_list)][summary_cols]
    cdr_sum = cdr_df.groupby('scen_id').sum()['2050']

    # Share of primary energy from renewables in 2030
    c1_prien_2030 = filtered_em.loc[
        filtered_em['Variable'] == 'Primary Energy']['2030']
    c1_renen_2030 = filtered_em.loc[
        filtered_em[
            'Variable'] == 'Primary Energy|Renewables (incl. Biomass)']['2030']
    c1_en_2030_df = pandas.DataFrame({
        'Primary Energy': c1_prien_2030,
        'Primary Energy|Renewables': c1_renen_2030})
    ren_share_2030 = (c1_en_2030_df['Primary Energy|Renewables'] /
        c1_en_2030_df['Primary Energy'])

    # Share of primary energy from renewables in 2050
    c1_prien_2050 = filtered_em.loc[
        filtered_em['Variable'] == 'Primary Energy']['2050']
    c1_renen_2050 = filtered_em.loc[
        filtered_em[
            'Variable'] == 'Primary Energy|Renewables (incl. Biomass)']['2050']
    c1_en_df = pandas.DataFrame({
        'Primary Energy': c1_prien_2050,
        'Primary Energy|Renewables': c1_renen_2050})
    ren_share_2050 = (c1_en_df['Primary Energy|Renewables'] /
        c1_en_df['Primary Energy'])

    # industrial process emissions
    ip_ccs = filtered_em.loc[
        filtered_em['Variable'] ==
        'Carbon Sequestration|CCS|Industrial Processes']['2050']

    filt_key_var = pandas.DataFrame({
        'AFOLU CO2e 2050': afolu_co2e,
        'Final energy 2050': fin_energy_2050,
        'Max yearly bioenergy': max_biom,
        'Cumulative CCS 2010-2050': cum_ccs,
        'CDR 2050': cdr_sum,
        'Ren share 2050': ren_share_2050,
        'Ren share 2030': ren_share_2030,
        'Industrial process CCS 2050': ip_ccs,
        'Cumulative gross fossil CO2 2020-2050': cum_gross_co2,
        })
    
    # add NZE key variables before calculating average across scenarios
    nze_key_var = pandas.read_csv(
        os.path.join(
            _PROJ_DIR, 'IEA/NZE_2023/NZE 2023 key variables summary.csv'))
    key_var_vals = pandas.concat([filt_key_var, nze_key_var])

    # average of scenarios in final scenario set
    key_var_avg = key_var_vals.mean()
    key_var_df = pandas.DataFrame({'scenario set mean': key_var_avg})
    key_var_df.to_csv(
        os.path.join(_OUT_DIR, 'scenario_set_key_var_summary.csv'))


def summarize_CO2():
    """Summarize gross EIP CO2 including harmonization with 2022 emissions."""
    # net EIP CO2 emissions, unharmonized
    ar6_key, ar6_scen = read_ar6_data()
    year_col = [col for col in ar6_scen if col.startswith('2')]
    c1_scen = ar6_key.loc[ar6_key['Category'] == 'C1']['scen_id']
    c1_filtered = sustainability_filters(c1_scen, ar6_scen)

    ar6_filled_em = fill_EIP_emissions(ar6_scen, c1_scen)
    net_co2_df = ar6_filled_em.loc[
        (ar6_filled_em['Variable'] ==
            'Emissions|CO2|Energy and Industrial Processes') &
        (ar6_filled_em['scen_id'].isin(c1_filtered))]
    net_co2_df.replace(0, numpy.nan, inplace=True)

    # gross EIP CO2 emissions (unharmonized)
    gross_em = calc_gross_eip_co2(ar6_filled_em, year_col)
    gross_co2_df = gross_em.loc[
        (gross_em['Variable'] ==
            'Emissions|CO2|Energy and Industrial Processes|Gross') &
        (gross_em['scen_id'].isin(c1_filtered))]
    gross_co2_df.replace(0, numpy.nan, inplace=True)

    # harmonize net CO2 emissions with historical
    hist_path = os.path.join(
        _PROJ_DIR, 'aneris_inputs', 'eip_co2_emissions_historical.csv')
    regions_path = os.path.join(
        _PROJ_DIR, 'aneris_inputs', 'regions_regions_sectors.csv')
    config_path = os.path.join(
        _PROJ_DIR, 'aneris_inputs', 'aneris_regions_sectors.yaml')
    aneris_dict = harmonize_to_historical(
        net_co2_df, hist_path, regions_path, config_path)

    # calculate gross EIP CO2 emissions from harmonized net emissions
    harm_df = aneris_dict['harmonized']
    harm_net_co2_df = harm_df.loc[
        harm_df['Variable'] == 'p|Emissions|CO2|Harmonized-DB']
    harm_net_co2_df['scen_id'] = harm_net_co2_df['Scenario']
    harm_net_co2_df['Variable'] = (
        'Emissions|CO2|Energy and Industrial Processes')
    non_net_df = ar6_filled_em.loc[
        (ar6_filled_em['Variable'] !=
            'Emissions|CO2|Energy and Industrial Processes') &
        (ar6_filled_em['scen_id'].isin(harm_df['Scenario'].unique()))]
    non_net_df.replace(0, numpy.nan, inplace=True)
    comb_df = pandas.concat([harm_net_co2_df, non_net_df])
    gross_harm_co2_df = calc_gross_eip_co2(comb_df, year_col, fill_df=True)
    gross_harm_co2_df['Variable'] = (
        'Emissions|CO2|Energy and Industrial Processes|Gross|Harmonized-DB')
    harm_net_co2_df['Variable'] = (
        'Emissions|CO2|Energy and Industrial Processes|Harmonized-DB')

    # gross EIP CO2 from IEA Net Zero Emissions by 2050 scenario
    iea_path = os.path.join(
        _PROJ_DIR, 'IEA', 'NZE_2023', 'gross_eip_co2_emissions.csv')
    iea_co2_df = pandas.read_csv(iea_path)
    iea_co2_df['scen_id'] = iea_co2_df['Source']
    iea_co2_df['Variable'] = (
        'Emissions|CO2|Energy and Industrial Processes|Gross')

    # add IEA to unharmonized and harmonized emissions df
    em_df = pandas.concat(
        [net_co2_df, harm_net_co2_df, gross_co2_df, gross_harm_co2_df,
        iea_co2_df])
    summary_cols = [
        str(idx) for idx in list(range(2020, 2051))] + ['Variable', 'scen_id']
    em_df = em_df[summary_cols]
    em_df.to_csv(os.path.join(_OUT_DIR, "combined_em_df_budg2050.csv"))


def cross_sector_benchmarks():
    """Calculate cross sector benchmarks from AR6 and IEA NZE scenarios."""
    # Median of filtered scenarios from AR6 database
    # This performs calculation of the scenario set for gross fossil
    # CO2 emissions, including harmonization with historical emissions within
    # a budget constraint calculated over 2022-2050.
    # The output is written to file at the path
    # os.path.join(_OUT_DIR, "combined_em_df_budg2050.csv")
    # summarize_CO2()

    # Gross fossil CO2 (output from above)
    harm_df = pandas.read_csv(
        os.path.join(_OUT_DIR, "combined_em_df_budg2050.csv"))

    # fill missing data via interpolation
    year_col = [col for col in harm_df if col.startswith('2')]
    harm_df.loc[
        harm_df['Variable'] ==
        'Emissions|CO2|Energy and Industrial Processes|Gross|Harmonized-DB',
        '2020'] = _2020_CO2
    harm_df['2021'] = numpy.nan
    gr_co2_df = harm_df.loc[
        (harm_df['Variable'] ==
        'Emissions|CO2|Energy and Industrial Processes|Gross|Harmonized-DB') |
        (harm_df['scen_id'] == 'NZE 2023')]
    gr_co2_df.replace(0, numpy.nan, inplace=True)
    gr_co2_df_interp = gr_co2_df[year_col].interpolate(axis=1)

    # median of harmonized scenarios, plus NZE
    gr_co2_med_ser = gr_co2_df_interp.quantile(q=0.5)
    gr_co2_med = pandas.DataFrame(gr_co2_med_ser).transpose()

    # summarize single-gas pathways for non-CO2 GHGs
    ar6_key, ar6_scen = read_ar6_data()
    year_col = [col for col in ar6_scen if col.startswith('2')]
    num_cols = ['2020', '2030', '2040', '2050']
    c1_scen = ar6_key.loc[ar6_key['Category'] == 'C1']['scen_id']

    eip_n2o_df = calc_eip_n2o(ar6_scen, year_col)
    eip_n2o_df.reset_index(inplace=True)
    med_c1_n2o = eip_n2o_df.loc[
        eip_n2o_df['scen_id'].isin(c1_scen)][num_cols].quantile(q=0.5)
    n2o_df = pandas.DataFrame(med_c1_n2o).transpose()
    n2o_df['Variable'] = 'Fossil N2O'
    n2o_df['Source'] = 'Cross sector benchmark'

    eip_ch4_df = calc_eip_ch4(ar6_scen, year_col)
    eip_ch4_df.reset_index(inplace=True)
    med_c1_ch4 = eip_ch4_df.loc[
        eip_ch4_df['scen_id'].isin(c1_scen)][num_cols].quantile(q=0.5)
    ch4_df = pandas.DataFrame(med_c1_ch4).transpose()
    ch4_df['Variable'] = 'Fossil CH4'
    ch4_df['Source'] = 'Cross sector benchmark'

    med_c1_hfc = ar6_scen.loc[
        (ar6_scen['Variable'] == 'Emissions|HFC') &
        (ar6_scen['scen_id'].isin(c1_scen)), ][num_cols].quantile(q=0.5)
    med_c1_pfc = ar6_scen.loc[
        (ar6_scen['Variable'] == 'Emissions|PFC') &
        (ar6_scen['scen_id'].isin(c1_scen)), ][num_cols].quantile(q=0.5)
    med_c1_sf6 = ar6_scen.loc[
        (ar6_scen['Variable'] == 'Emissions|SF6') &
        (ar6_scen['scen_id'].isin(c1_scen)), ][num_cols].quantile(q=0.5)

    # add CO2eq from non-CO2 GHGs
    non_co2_df = pandas.DataFrame(
        {'N2O (Mt CO2e)': med_c1_n2o * _N2O_GWP100_AR6 * _KT_to_MT,
        'CH4 (Mt CO2e)': med_c1_ch4 * _CH4FOSS_GWP100_AR6,
        'HFC (Mt CO2e)': med_c1_hfc * _HFC_GWP100_AR6 * _KT_to_MT,
        'PFC (Mt CO2e)': med_c1_pfc * _PFC_GWP100_AR6 * _KT_to_MT,
        'SF6 (Mt CO2e)': med_c1_sf6 * _SF6_GWP100_AR6 * _KT_to_MT}).transpose()
    non_co2_df['Variable'] = non_co2_df.index
    non_co2_df['Source'] = 'Cross sector benchmark'

    co2_df = gr_co2_med[num_cols]
    co2e_df = pandas.concat([co2_df, non_co2_df])
    sum_ser = co2e_df[num_cols].sum()
    sum_df = pandas.DataFrame(sum_ser).transpose()
    sum_df['Variable'] = 'Gross fossil CO2e'
    sum_df['Source'] = 'Cross sector benchmark'

    # summary for inclusion in the report
    co2_df['Variable'] = 'Gross fossil CO2'
    co2_df['Source'] = 'Cross sector benchmark'
    summary_df = pandas.concat([co2_df, sum_df, non_co2_df])
    summary_df['percch_2030'] = (
        summary_df['2030'] - summary_df['2020']) / summary_df['2020']
    summary_df['percch_2040'] = (
        summary_df['2040'] - summary_df['2020']) / summary_df['2020']
    summary_df['percch_2050'] = (
        summary_df['2050'] - summary_df['2020']) / summary_df['2020']
    summary_df.to_csv(
        os.path.join(_OUT_DIR, '20240511_cs_benchmark_summary.csv'),
        index=False)


def main():
    summarize_c1_key_var()
    summarize_filtered_key_var()
    cross_sector_benchmarks()


if __name__ == '__main__':
    main()