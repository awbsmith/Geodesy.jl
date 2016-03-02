### transformations for dynamic datums.  Can Proj4 do this?

#=TODO: make this
function convert{T,Y}(::Type{StaticLLA{GDA94}}, X::StaticLLA{T, Y})


    // ITRF transformation parameters from
    // J. Dawson and A. Woods, "ITRF to GDA94 coordinate transforms",
    // Journal of Applied Geodesy, 4, p. 189 (2010):
    //
    // Format of the following table is:
    //   tx ty tz sc rx ry rz (value at 1994 reference epoch)
    //   tx ty tz sc rx ry rz (derivatives)
    //
    // Units are:
    //   mm    mm    mm    ppb    mas    mas    mas
    //   mm/yr mm/yr mm/yr ppb/yr mas/yr mas/yr mas/yr
    //
    // where ppb = parts per billion; mas = milli arc seconds; yr = year
    //

    const double* params = 0;
    const double* rates = 0;
    const double itrf2008_params[7] = {-84.68, -19.42, 32.01, 9.710, -0.4254, 2.2578, 2.4015};
    const double itrf2008_rates [7] = {1.42, 1.34, 0.90, 0.109, 1.5461, 1.1820, 1.1551};
    const double itrf2005_params[7] = {-79.73, -6.86, 38.03, 6.636, -0.0351, 2.1211, 2.1411};
    const double itrf2005_rates [7] = {2.25, -0.62, -0.56, 0.294, 1.4707, 1.1443, 1.1701};
    const double itrf2000_params[7] = {-45.91, -29.85, -20.37, 7.070, -1.6705, 0.4594, 1.9356};
    const double itrf2000_rates [7] = {-4.66, 3.55, 11.24, 0.249, 1.7454, 1.4868, 1.2240};
    const double itrf1997_params[7] = {-14.63, -27.62, -25.32, 6.695, -1.7893, -0.6047, 0.9962};
    const double itrf1997_rates [7] = {-8.60, 0.36, 11.25, 0.007, 1.6394, 1.5198, 1.3801};
    const double itrf1996_params[7] = {24.54, -36.43, -68.12, 6.901, -2.7359, -2.0431, 0.3731};
    const double itrf1996_rates [7] = {-21.80, 4.71, 26.27, 0.388, 2.0203, 2.1735, 1.6290};
    if      (inputDatumName == "ITRF2008") { params = itrf2008_params; rates = itrf2008_rates; }
    else if (inputDatumName == "ITRF2005") { params = itrf2005_params; rates = itrf2005_rates; }
    else if (inputDatumName == "ITRF2000") { params = itrf2000_params; rates = itrf2000_rates; }
    else if (inputDatumName == "ITRF1997") { params = itrf1997_params; rates = itrf1997_rates; }
    else if (inputDatumName == "ITRF1996") { params = itrf1996_params; rates = itrf1996_rates; }
    else
    {
        throw RoamesError("Transform from input datum %s not supported", inputDatum);
    }
    if (!inputDatumHadTime)
    {
        throw RoamesError("Date must be provided for %s, "
                          "since the datum is time dependent", inputDatum);
    }

    // Unit conversion factors to get into meters & radians
    const double mas2rad = 1e-3/(60*60) * M_PI/180;
    const double unitConv[7] = {1e-3, 1e-3, 1e-3, 1e-9, mas2rad, mas2rad, mas2rad};
    // time is in fractional years, relative to 1994.
    double t = inputDatumTime - 1994;
    double coeffs[7] = {0};
    for (int i = 0; i < 7; ++i)
        coeffs[i] = unitConv[i] * (params[i] + rates[i]*t);
    m_Tx = coeffs[0];
    m_Ty = coeffs[1];
    m_Tz = coeffs[2];
    m_Sc = coeffs[3];
    m_Rx = coeffs[4];
    m_Ry = coeffs[5];
    m_Rz = coeffs[6];
    LOG_DEBUG("%s -> %s tansformation parameters: "
              "T = (%.5f %.5f %.5f), Sc = %.5e, R = (%.5e, %.5e, %.5e)",
              inputDatum, outputDatum, m_Tx, m_Ty, m_Tz, m_Sc, m_Rx, m_Ry, m_Rz);

=#
