# Analysis-of-in-vivo-skin-properties-using-elastic-wave-measurements
Github to Accompany the publication of the same name


The datasets containing the averaged Reviscometer data for the natural and stretched configurations were titled:
"Reviscometer Avg Data - RRT Units.csv"
"Reviscometer Avg Data - Stretched - RRT Units.csv"

These csv files had the structure:


| Subject_Id | Gender | Age | 0 Degrees | 10 Degrees | 20 Degrees | ... | 350 Degrees | Min_RRT | Max_RRT | Range_RRT | Max_Min_Ratio |
|------------|--------|-----|-----------|------------|------------|-----|-------------|---------|---------|-----------|---------------|
| 1          | Female | 3   | -         | -          | -          | ... | -           | -       | -       | -         | -             |
| 2          | Male   | 3   | -         | -          | -          | ... | -           | -       | -       | -         | -             |
| 3          | Female | 5   | -         | -          | -          | ... | -           | -       | -       | -         | -             |



Ellipses were fit to these datasets and the parameters were extracted using the code in:
"Fit Ellipses to Reviscometer Data"

The resulting csv (used for modelling) had the structure:

|     | Subject_ID | Gender | Age | Centre_x | Centre_y | Axis_1 | Axis_2 | Semi-major_Axis | Semi-minor_Axis | Angle_Rads | Angle_Deg | Angle_to_semi-maj_Rads | Angle_to_semi-maj_Deg | Eccentricity | RSS | RSS_Area_Norm | Area | Configuration |
|-----|------------|--------|-----|----------|----------|--------|--------|-----------------|-----------------|------------|-----------|------------------------|-----------------------|--------------|-----|---------------|------|---------------|
| 1   | 1          | Female | 3   | -        | -        | -      | -      | -               | -               | -          | -         | -                      | -                     | -            | -   | -             | -    | Natural       |
| 2   | 2          | Male   | 3   | -        | -        | -      | -      | -               | -               | -          | -         | -                      | -                     | -            | -   | -             | -    | Natural       |
| 3   | 3          | Female | 5   | -        | -        | -      | -      | -               | -               | -          | -         | -                      | -                     | -            | -   | -             | -    | Natural       |
| ... | ...        | ...    | ... | ...      | ...      | ...    | ...    | ...             | ...             | ...        | ...       | ...                    | ...                   | ...          | ... | ...           | ...  | ...           |
| 1   | 1          | Female | 3   | -        | -        | -      | -      | -               | -               | -          | -         | -                      | -                     | -            | -   | -             | -    | Stretched     |
| 2   | 2          | Male   | 3   | -        | -        | -      | -      | -               | -               | -          | -         | -                      | -                     | -            | -   | -             | -    | Stretched     |
| 3   | 3          | Female | 5   | -        | -        | -      | -      | -               | -               | -          | -         | -                      | -                     | -            | -   | -             | -    | Stretched     |
