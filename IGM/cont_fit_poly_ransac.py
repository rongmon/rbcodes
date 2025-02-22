import numpy as np
from sklearn.linear_model import RANSACRegressor, LinearRegression
from sklearn.preprocessing import PolynomialFeatures
import matplotlib.pyplot as plt

def fit_polynomial_ransac(wave, flux, error, degree, residual_threshold=0.5, n_bootstrap=100):
    """
    Fits a polynomial using RANSAC to remove outliers and computes the uncertainty of the fit using bootstrap resampling.
    
    Parameters:
        wave (array): Wavelength values.
        flux (array): Observed flux values.
        error (array): Uncertainty in flux values.
        degree (int): Degree of the polynomial to fit.
        residual_threshold (float): Threshold for identifying outliers based on residuals.
        n_bootstrap (int): Number of bootstrap resampling iterations.
    
    Returns:
        flux_fit (array): Best-fit polynomial evaluated at wave points.
        model_error (array): Estimated uncertainty in the fit (standard deviation of bootstrap predictions).
    """
    
    # Transform the wave values into polynomial features
    poly = PolynomialFeatures(degree)
    wave_poly = poly.fit_transform(wave.reshape(-1, 1))
    
    # Apply RANSAC with a Linear Regression estimator (since RANSACRegressor needs a regression model)
    ransac = RANSACRegressor(estimator=LinearRegression(), residual_threshold=residual_threshold)
    ransac.fit(wave_poly, flux)

    # Extract inliers
    inliers = ransac.inlier_mask_
    X_inliers, y_inliers, err_inliers = wave[inliers].reshape(-1, 1), flux[inliers], error[inliers]
    
    # Fit the model using inliers only (weighted by error)
    linreg = LinearRegression()
    linreg.fit(poly.transform(X_inliers), y_inliers)
    
    # Get the fitted flux for all data points from the initial fit
    flux_fit = linreg.predict(poly.transform(wave.reshape(-1, 1)))

    # Bootstrap resampling to estimate model uncertainty
    bootstrap_predictions = np.zeros((n_bootstrap, len(wave)))
    for i in range(n_bootstrap):
        # Resample with replacement
        indices = np.random.choice(len(wave), size=len(wave), replace=True)
        X_bootstrap = wave[indices].reshape(-1, 1)
        y_bootstrap = flux[indices]
        
        # Fit the model on the bootstrap sample
        linreg.fit(poly.transform(X_bootstrap), y_bootstrap)
        
        # Predict on the original data points
        bootstrap_predictions[i, :] = linreg.predict(poly.transform(wave.reshape(-1, 1)))

    # Compute the standard deviation of the bootstrap predictions (model uncertainty)
    model_error = np.std(bootstrap_predictions, axis=0)

    return flux_fit, model_error
