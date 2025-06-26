# correlation_power_analysis

A Python module for conducting power analysis in correlational studies. This tool helps researchers determine appropriate sample sizes, calculate statistical power, and find minimum detectable effect sizes for Pearson correlation analyses.

## Installation

Simply download `correlation_power.py` and place it in your project directory, or clone this repository:

```bash
git clone https://github.com/yourusername/correlation-power.git
```

## Dependencies

```bash
pip install numpy scipy matplotlib pandas
```

## Quick Start

```python
from correlation_power import CorrelationPowerAnalysis, quick_sample_size, quick_power

# Quick calculations
n_needed = quick_sample_size(r=0.3, power=0.8)  # Sample size for r=0.3, 80% power
power = quick_power(n=100, r=0.25)  # Power with N=100, r=0.25

# Comprehensive analysis
analyzer = CorrelationPowerAnalysis()
results = analyzer.analyze_study(
    target_r=0.3,           # Expected correlation
    num_tests=3,            # Number of correlations (for multiple comparisons)
    power=0.8,              # Desired power
    show_plot=True          # Show power curves
)
```

## Main Functions

### Sample Size Calculation
```python
analyzer = CorrelationPowerAnalysis()
n = analyzer.sample_size(r=0.3, power=0.8, alpha=0.05)
print(f"Required sample size: {n}")
```

### Power Calculation
```python
power = analyzer.power(n=100, r=0.25, alpha=0.05)
print(f"Statistical power: {power:.3f}")
```

### Minimum Detectable Effect
```python
min_r = analyzer.minimum_effect(n=80, power=0.8)
print(f"Minimum detectable effect: {min_r:.3f}")
```

### Comprehensive Study Analysis
```python
results = analyzer.analyze_study(
    target_r=0.3,     # Expected correlation
    num_tests=3,      # Number of tests for multiple comparisons
    power=0.8,        # Desired power
    show_plot=True    # Show power curve plots
)
```

## Visualizations

### Power Curves
```python
analyzer.plot_power_curves(
    effect_sizes=[0.2, 0.3, 0.4, 0.5],
    n_range=(20, 150),
    target_power=0.8
)
```

### Power Table
```python
table = analyzer.power_table(
    effect_sizes=[0.2, 0.3, 0.4, 0.5],
    sample_sizes=[50, 75, 100, 125, 150]
)
print(table)
```

## Examples

### Example 1: Basic Power Analysis
```python
from correlation_power import CorrelationPowerAnalysis

analyzer = CorrelationPowerAnalysis()

# How many participants do I need for r=0.3 with 80% power?
n = analyzer.sample_size(r=0.3, power=0.8)
print(f"Required sample size: {n}")

# What power do I have with 100 participants for r=0.25?
power = analyzer.power(n=100, r=0.25)
print(f"Statistical power: {power:.3f}")
```

### Example 2: Multiple Comparisons
```python
# Study with 3 correlations
results = analyzer.analyze_study(
    target_r=0.3,
    num_tests=3,  # Bonferroni correction applied
    power=0.8
)
```

### Example 3: Compare Different Studies
```python
from correlation_power import compare_studies

studies = [
    {'name': 'Study A', 'r': 0.3, 'n': 85, 'num_tests': 1},
    {'name': 'Study B', 'r': 0.25, 'n': 100, 'num_tests': 3},
    {'name': 'Study C', 'r': 0.4, 'n': 60, 'num_tests': 2}
]

comparison = compare_studies(studies)
```

## Quick Functions

For simple calculations without creating an analyzer object:

```python
from correlation_power import quick_sample_size, quick_power, quick_min_effect

# Sample size needed
n = quick_sample_size(r=0.3, power=0.8)

# Power achieved  
power = quick_power(n=100, r=0.25)

# Minimum detectable effect
min_r = quick_min_effect(n=80, power=0.8)
```

## Features

- ✅ Sample size calculation for given effect size and power
- ✅ Power calculation for given sample size and effect size
- ✅ Minimum detectable effect size calculation
- ✅ Multiple comparisons correction (Bonferroni)
- ✅ Power curve visualization
- ✅ Comprehensive study analysis with detailed reporting
- ✅ Quick utility functions for common calculations
- ✅ Input validation and error handling
- ✅ Support for one-tailed and two-tailed tests

## Method Details

### CorrelationPowerAnalysis Class

| Method | Description | Parameters |
|--------|-------------|------------|
| `sample_size()` | Calculate required N | `r`, `power`, `alpha`, `two_tailed` |
| `power()` | Calculate statistical power | `n`, `r`, `alpha`, `two_tailed` |
| `minimum_effect()` | Find min detectable correlation | `n`, `power`, `alpha`, `two_tailed` |
| `analyze_study()` | Comprehensive analysis | `target_r`, `n`, `power`, `alpha`, `num_tests`, `show_plot` |
| `plot_power_curves()` | Visualize power curves | `effect_sizes`, `n_range`, `alpha`, `target_power` |
| `power_table()` | Create power table | `effect_sizes`, `sample_sizes`, `alpha` |

### Quick Functions

| Function | Description | Parameters |
|----------|-------------|------------|
| `quick_sample_size()` | Fast sample size calculation | `r`, `power`, `alpha`, `two_tailed` |
| `quick_power()` | Fast power calculation | `n`, `r`, `alpha`, `two_tailed` |
| `quick_min_effect()` | Fast min effect calculation | `n`, `power`, `alpha`, `two_tailed` |
| `compare_studies()` | Compare multiple studies | `studies`, `show_table` |

## Citation

If you use this tool in your research, please cite:

```
[Your Name] (2024). Correlation Power Analysis: A Python module for power analysis in correlational studies. 
GitHub: https://github.com/yourusername/correlation-power
```

## License

MIT License - feel free to use, modify, and distribute.

## Contributing

Contributions welcome! Please feel free to submit issues or pull requests.
