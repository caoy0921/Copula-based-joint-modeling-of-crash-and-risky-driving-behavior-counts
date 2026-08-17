# Copula-based-joint-modeling-of-crash-and-risky-driving-behavior-counts
R code for jointly modeling crash counts and risky driving behavior counts using Copula models for freeway safety evaluation.
The complete analysis is provided in: copula_analysis.R
## Data availability
The original crash and risky driving behavior data used in this study cannot be publicly shared because they are subject to confidentiality and third-party data restrictions. Therefore, the original dataset is not included in this repository. Researchers who wish to apply the code should prepare their own input data according to the variable definitions provided below and save the file as `data/example_data.xlsx`.
The input dataset should contain the following variables:
- `TIME`: date and hour of observation;
- `road`: freeway segment;
- `FC`: crash count;
- `FB`: risky driving behavior count;
- `V`: traffic volumn;
- `S`: traffic speed;
- `LH`: segment length;
- `S_std`: standard deviation of speed;
- `EX`: ramp;
- `LE`: road alignment.
