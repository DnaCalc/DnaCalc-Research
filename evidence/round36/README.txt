Round 36: original scalar questions and bounded arithmetic observations

All retained answers are from Excel 16.0 build 20527 through the question
service. files.json binds each compressed answer and its producer. Producers
use the paths from the original run; change those paths to replay elsewhere.
The first cumulative timing request included a subnormal and negative zero;
the service rejected its input audit, and it contributes no observations.

LINEST one-predictor and exactly collinear two-predictor fits use symmetric
responses, shifts, predictor scales, reversed rows, and both const values.
They contain 480 scalar selections. The intercept without a constant is an
error in all 48 relevant selections. On 24 eligible intercept fits per design,
the ordinary centred sum-of-squares expression matches 22 and fails two;
the two failed one-predictor cases and their exactly collinear counterparts
return 41143d1362487ba7 versus the expression's 41143d1362487ba8.
This compares complete published values and does not expose a private QR leaf.

An independent adjacent-input run takes 600 selections from 60 LINEST and
60 LOGEST fits. The shifts straddle binary64 values near 1, 16 and 1024;
both row orders and two predictor scales are present. Centred residual-sum
expressions match 44/60 LINEST and 41/60 LOGEST; expressions using the
published slope standard error match 36/60 and 40/60. Neither public result
identifies the hidden inverse-R or residual arithmetic. These are refutations
of the declared expressions, not an exhaustive regression search.

On coupon dates with zero yield, 150 DURATION and 150 MDURATION outputs
match each other and each of 16 simple strict summation labels on the initial
dyadic controls. They do not distinguish their sum order. Before asking further,
108 additional cases were chosen where eight fixed labels predict at least
three different values. Only the forward coupon sum, weighted coupon terms,
and redemption added last matches all 108. A separate 108-case holdout with
different periods and coupon exponents matches that expression for both
functions. With yield equal to frequency, every discount is an exact power of
two; another 108 cases match the same forward sum for each function. These
are on-coupon controls only; fractional and general discount stages remain open.

For the 25 original long cumulative ranges, fresh whole and two-part readings
reproduce all 50 previous whole-range answers and are additive after a
binary64 addition on 10 CUMPRINC and 11 CUMIPMT cases. The remaining cases
do not identify the internal producer or summation order.
Sixteen payment-timing inputs per function show only exact zero and one
accepted, including rejection of 1e-13 and the two neighbours of one. The
qualified input readbacks and original output errors are retained.

For 120 POISSON.DIST CDF controls around four exponent or dispatch positions,
the paired CHISQ.DIST.RT outputs agree on all 120; subtracting the paired
GAMMA.DIST CDF from one agrees on 76. This establishes the observed public
relation on these inputs only; the internal Poisson/gamma stages remain open.
The prior six number-to-text failures had no new separating question, so no
additional conversion model was selected. The problem list's INDEX examples
include observers of SORT spills and formatted HSTACK results, so the label
alone is not evidence of an INDEX arithmetic defect.
