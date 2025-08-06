*Overall, the manuscript provides a solution for an important task – averaging of ESM ensembles -
that is both well-grounded in theory and easy to apply. The (sometimes strong) assumptions are
clearly stated, and the method overall is well presented in terms of derivation and examples.*

*A few questions remain:
1) Calculation of models weights. In eqs (7) – (15), I take it that index “i” goes over all models
except the best performing model (“m”), which serves as the reference. Does that mean the
best performing model will not be included in the weighted average? Please clarify.*

The index *i* runs over all models, including the best performing model *m*. 
For the calculation of the weights, the best performing model has a value $\Delta_m =0$, which
implies that the value of the numerator in the computation of the weights (eq. 13) is equal to 1,
and for all other models the value in the numerator is less than 1. 

*2) Patterns of model performance in space and time (see also line 180): The authors resolve
spatial patterns of model performance by calculating model weights grid-by-grid over all
points in time. This inherently assumes temporal invariance of relative model performance,
which clearly is not the case. Therefore I wonder if it would not be more appropriate to
derive model residuals (and deltas) from space-time regions rather than time-regions alone.
Please comment.*

In principle we agree in that model performance should account for spatial and temporal 
covariations. However, it is not trivial to include these covariations in our information-
theoretic approach because they have to be treated in the context of mutual information. 
The same reasoning applies for covariation among different models, which should be 
treated as mutual information and not just quantifying a covariance matrix. 
We believe that this is a topic that deserves further investigation, and we added a paragraph
in the Discussion section addressing this topic. 

*3) In Sect. 4, simple averaging is used as a benchmark to compare the weighted averaging
proposed by the authors. While e.g. Fig. 5 clearly show differences among the methods, it is
not clear if the weighted average really provides the better (in terms of smaller disagreement
from observations) estimate that the simple average. I am quite sure this will be the case,
but please add and discuss the related numbers.*

From the theoretical point of view, the inverse-variance weighted average is the most efficient
estimator of the mean under the Maximum Likelihood estimation theory. The numerical results
support this claim showing a much reduced variance in comparison to the weighted average, and
narrower prediction intervals.  

*My overall recommendation is to publish after these minor points have been suitably addressed.*



