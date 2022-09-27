# Circular-design-for-total-effects
This project creates optimal circular designs for estimaing total effects as supplementary materials for "Circular designs for total effects under interference models"

Many thanks for your careful reading. We are now uploading the R precedures for generating optimal designs in this paper. They are separated into 4 different files. The sequence.R file generates equivalence clasess; the information.R file calculates the information matrix of equivalence classes; the design.R file generates optimal/efficient exact designs; the efficiency.R file gives the efficiency of the exact designs generated above. This 4 files should be compiled sequentially and the main function is find good design(). Please, refer to the R files for more details.

It should be mentioned that running this code needs a special R package gurobi, which needs registration and installation from gurobi IBM.
