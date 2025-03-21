import numpy as np
from sklearn.linear_model import RANSACRegressor, LinearRegression
from sklearn.preprocessing import PolynomialFeatures

def fit_polynomial_ransac(wave, flux, error, degree, residual_threshold=0.5, n_bootstrap=False, verbose=True, return_model=False):
    """
    Fits a polynomial using RANSAC to remove outliers and estimates uncertainty either by bootstrap resampling
    or by computing standard error if bootstrap is disabled.
    
    Parameters:
        wave (array): Wavelength values.
        flux (array): Observed flux values.
        error (array): Uncertainty in flux values.
        degree (int): Degree of the polynomial to fit.
        residual_threshold (float): Threshold for identifying outliers based on residuals.
        n_bootstrap (int or bool): Number of bootstrap resampling iterations, or False to skip bootstrap.
        verbose (bool): If True, prints warnings about poor fit or insufficient inliers.
        return_model (bool): If True, returns a function to evaluate the fitted model at new wavelengths.
    
    Returns:
        flux_fit (array): Best-fit polynomial evaluated at wave points.
        model_error (array): Estimated uncertainty in the fit.
        model (callable, optional): Function to compute the continuum model at new wavelength values.
    """
    # Transform the wave values into polynomial features
    poly = PolynomialFeatures(degree)
    wave_poly = poly.fit_transform(wave.reshape(-1, 1))
    np.random.seed(42)

    # Apply RANSAC with a Linear Regression estimator
    ransac = RANSACRegressor(base_estimator=LinearRegression(), residual_threshold=residual_threshold, random_state=42)
    ransac.fit(wave_poly, flux)

    # Extract inliers
    inliers = ransac.inlier_mask_
    num_inliers = np.sum(inliers)
    inlier_ratio = num_inliers / len(wave)

    if verbose:
        if num_inliers == 0:
            print("Warning: RANSAC failed to find any inliers. Fit is unreliable.")
        elif inlier_ratio < 0.3:
            print(f"Warning: Only {inlier_ratio:.2%} of data points classified as inliers. Fit may be poor.")

    # If no inliers are found, return NaN arrays
    if num_inliers == 0:
        return np.full_like(wave, np.nan), np.full_like(wave, np.nan), (lambda x: np.full_like(x, np.nan)) if return_model else None

    X_inliers, y_inliers, err_inliers = wave[inliers].reshape(-1, 1), flux[inliers], error[inliers]
    
    # Fit the model using inliers only
    linreg = LinearRegression()
    linreg.fit(poly.transform(X_inliers), y_inliers)
    
    # Get the fitted flux for all data points
    flux_fit = linreg.predict(poly.transform(wave.reshape(-1, 1)))
    
    if n_bootstrap:
        # Bootstrap resampling to estimate model uncertainty
        bootstrap_predictions = np.zeros((n_bootstrap, len(wave)))
        for i in range(n_bootstrap):
            indices = np.random.choice(len(wave), size=len(wave), replace=True)
            X_bootstrap, y_bootstrap = wave[indices].reshape(-1, 1), flux[indices]
            
            # Fit the model on the bootstrap sample
            linreg.fit(poly.transform(X_bootstrap), y_bootstrap)
            
            # Predict on the original data points
            bootstrap_predictions[i, :] = linreg.predict(poly.transform(wave.reshape(-1, 1)))
        
        # Compute the standard deviation of bootstrap predictions (uncertainty)
        model_error = np.std(bootstrap_predictions, axis=0)
    
    else:
        # Compute standard error of the fit using inliers' residuals
        residuals = y_inliers - linreg.predict(poly.transform(X_inliers))
        
        # Standard error per point using error propagation
        model_error = np.sqrt(np.mean(residuals**2)) * np.ones_like(flux_fit) / np.sqrt(num_inliers)
    
    # Define the model function for new wavelength values
    def model(wave_new):
        return linreg.predict(poly.transform(np.array(wave_new).reshape(-1, 1)))
    
    if return_model:
        return flux_fit, model_error, model
    else:
        return flux_fit, model_error
