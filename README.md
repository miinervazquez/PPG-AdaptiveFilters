# PPG-AdaptiveFilters
Code to reproduce the results of the article:  
**"Spectral coherence as a method for evaluating adaptive filtering techniques in photoplethysmography"**

## Manuscript Information
- **Manuscript Title:** Spectral coherence as a method for evaluating adaptive filtering techniques in photoplethysmography  
- **Submission ID:** 9211  
- **Authors:** [Minerva Vázquez], [J. Gerardo Ávalos Ochoa], [J. Carlos Sánchez García], [Brayans Becerra Luna].
  
## Description
This repository contains MATLAB code for calculating spectral coherence and power spectral density, as well as implementing adaptive filters used in the study.

### Code Structure
1. **Main filter implementation files:**
   - `filtrosLMS.m`: Contains the implementation of the LMS family of filters.
   - `filtrosAPA.m`: Contains the implementation of the APA family of filters.

   These files call auxiliary functions (`.m` files) that perform specific calculations for each filter. The names of these functions correspond to the abbreviation of the algorithm they implement, making them easy to identify and use.

2. **Experimental data:**
   - The data from the seven volunteers is stored in `.xlsx` files.
   - Each file contains 18 columns with the following structure:
     1. Sample number
     2. PPG signal with 60 Hz noise
     3. PPG signal filtered to remove 60 Hz noise
     4. Contaminated PPG signal
     5. Motion signal in the X-axis
     6. Motion signal in the Y-axis
     7. Motion signal in the Z-axis
     8. Empty column
     9. Signal filtered with the LMS filter
     10. Signal filtered with the SLMS filter
     11. Empty column
     12. Signal filtered with the SSLMS filter
     13. Signal filtered with the NLMS filter
     14. Signal filtered with the VSSNLMS filter
     15. Empty column
     16. Signal filtered with the APA filter
     17. Signal filtered with the VSSAPA filter
     18. Signal filtered with the VAPA filter
