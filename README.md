# PPLS
Please also install OmicsPLS (from selbouhaddani/OmicsPLS) and ggplot2 in order to make PPLS work.

```
require(devtools)
install_github("selbouhaddani/PPLS/Package/PPLS")
library(PPLS)
```

If this doesn't work, try downloading the tar or zip file and install from local disk.

For questions mail me or file an issue

### Example
Simultaneous PPLS algorithm: `PPLS_simult(X, Y, nr_comp, ...)`
