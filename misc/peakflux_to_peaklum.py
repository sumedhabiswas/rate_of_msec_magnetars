import numpy as np
from astropy.cosmology import Planck18 as cosmo 

def calculate_luminosity(flux, redshift):
    # Convert redshift to luminosity distance in cm
    distance_cm = cosmo.luminosity_distance(redshift).to('cm').value
    
    # Calculate luminosity
    luminosity = 4 * np.pi * distance_cm**2 * flux
    return luminosity


if __name__ == "__main__":
    # Get the flux from the user
    while True:
        try:
            flux = float(input("Enter the peak unabsorbed flux in erg/s/cm^2: "))
            break
        except ValueError:
            print("Invalid input. Please enter a numeric value for the flux.")

    # Redshift range
    redshifts = np.linspace(0.5, 1.5, 10)  # 10 points between z=0.5 and z=1.5

    # Calculate luminosities for each redshift
    luminosities = [calculate_luminosity(flux, z) for z in redshifts]


    for z, lum in zip(redshifts, luminosities):
        print(f"Redshift: {z:.2f}, Peak Luminosity: {lum:.2e} erg/s")


    average_luminosity = np.mean(luminosities)
    print(f"\nAverage Peak Luminosity across redshifts 0.5 to 1.5: {average_luminosity:.2e} erg/s")

