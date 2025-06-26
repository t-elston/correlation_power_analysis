"""
Correlation Power Analysis Tool
A general-purpose tool for power analysis in correlational studies

Author: Thomas Elston (elston@utexas.edu)
License: MIT
"""

import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import pandas as pd
from typing import Union, List, Optional

class CorrelationPowerAnalysis:
    """
    A class for conducting power analysis for correlation studies.
    
    This tool can calculate:
    - Required sample size for a given effect size and power
    - Statistical power for a given sample size and effect size
    - Minimum detectable effect size for a given sample size and power
    """
    
    def __init__(self):
        pass
    
    @staticmethod
    def fisher_z_transform(r: float) -> float:
        """Convert correlation r to Fisher's Z"""
        return 0.5 * np.log((1 + r) / (1 - r))
    
    @staticmethod
    def inverse_fisher_z(z: float) -> float:
        """Convert Fisher's Z back to correlation"""
        return (np.exp(2 * z) - 1) / (np.exp(2 * z) + 1)
    
    def power(self, n: int, r: float, alpha: float = 0.05, two_tailed: bool = True) -> float:
        """
        Calculate statistical power for correlation test
        
        Parameters:
        -----------
        n : int
            Sample size
        r : float
            Expected correlation coefficient
        alpha : float, default 0.05
            Significance level
        two_tailed : bool, default True
            Whether to use two-tailed test
            
        Returns:
        --------
        float
            Statistical power (0-1)
        """
        z_r = self.fisher_z_transform(abs(r))
        se = 1 / np.sqrt(n - 3)
        
        if two_tailed:
            z_critical = stats.norm.ppf(1 - alpha/2)
        else:
            z_critical = stats.norm.ppf(1 - alpha)
            
        z_score = z_r / se
        
        if two_tailed:
            power = 1 - stats.norm.cdf(z_critical - z_score) + stats.norm.cdf(-z_critical - z_score)
        else:
            power = 1 - stats.norm.cdf(z_critical - z_score)
            
        return min(power, 1.0)  # Cap at 1.0
    
    def sample_size(self, r: float, power: float = 0.8, alpha: float = 0.05, 
                   two_tailed: bool = True) -> int:
        """
        Calculate required sample size for correlation test
        
        Parameters:
        -----------
        r : float
            Expected correlation coefficient
        power : float, default 0.8
            Desired statistical power
        alpha : float, default 0.05
            Significance level
        two_tailed : bool, default True
            Whether to use two-tailed test
            
        Returns:
        --------
        int
            Required sample size
        """
        z_r = self.fisher_z_transform(abs(r))
        
        if two_tailed:
            z_alpha = stats.norm.ppf(1 - alpha/2)
        else:
            z_alpha = stats.norm.ppf(1 - alpha)
            
        z_beta = stats.norm.ppf(power)
        
        n = ((z_alpha + z_beta) / z_r) ** 2 + 3
        
        return int(np.ceil(n))
    
    def minimum_effect(self, n: int, power: float = 0.8, alpha: float = 0.05, 
                      two_tailed: bool = True) -> float:
        """
        Calculate minimum detectable effect size
        
        Parameters:
        -----------
        n : int
            Sample size
        power : float, default 0.8
            Desired statistical power
        alpha : float, default 0.05
            Significance level
        two_tailed : bool, default True
            Whether to use two-tailed test
            
        Returns:
        --------
        float
            Minimum detectable correlation
        """
        if two_tailed:
            z_alpha = stats.norm.ppf(1 - alpha/2)
        else:
            z_alpha = stats.norm.ppf(1 - alpha)
            
        z_beta = stats.norm.ppf(power)
        se = 1 / np.sqrt(n - 3)
        
        z_required = (z_alpha + z_beta) * se
        min_r = self.inverse_fisher_z(z_required)
        
        return min_r
    
    def power_table(self, effect_sizes: List[float], sample_sizes: List[int], 
                   alpha: float = 0.05, two_tailed: bool = True) -> pd.DataFrame:
        """
        Create a power table for multiple effect sizes and sample sizes
        
        Parameters:
        -----------
        effect_sizes : List[float]
            List of correlation coefficients to test
        sample_sizes : List[int]
            List of sample sizes to test
        alpha : float, default 0.05
            Significance level
        two_tailed : bool, default True
            Whether to use two-tailed test
            
        Returns:
        --------
        pd.DataFrame
            Power table with effect sizes as columns and sample sizes as rows
        """
        power_matrix = []
        
        for n in sample_sizes:
            power_row = []
            for r in effect_sizes:
                power_val = self.power(n, r, alpha, two_tailed)
                power_row.append(power_val)
            power_matrix.append(power_row)
        
        df = pd.DataFrame(power_matrix, 
                         index=[f'N={n}' for n in sample_sizes],
                         columns=[f'r={r:.2f}' for r in effect_sizes])
        
        return df.round(3)
    
    def plot_power_curves(self, effect_sizes: List[float], n_range: tuple = (10, 200), 
                         alpha: float = 0.05, target_power: float = 0.8,
                         figsize: tuple = (10, 6)):
        """
        Plot power curves for different effect sizes
        
        Parameters:
        -----------
        effect_sizes : List[float]
            List of correlation coefficients to plot
        n_range : tuple, default (10, 200)
            Range of sample sizes to plot (min, max)
        alpha : float, default 0.05
            Significance level
        target_power : float, default 0.8
            Target power line to draw
        figsize : tuple, default (10, 6)
            Figure size
        """
        sample_sizes = np.arange(n_range[0], n_range[1] + 1, 2)
        
        plt.figure(figsize=figsize)
        
        for r in effect_sizes:
            power_values = [self.power(n, r, alpha) for n in sample_sizes]
            plt.plot(sample_sizes, power_values, label=f'r = {r:.2f}', linewidth=2)
        
        plt.axhline(y=target_power, color='red', linestyle='--', alpha=0.7, 
                   label=f'{int(target_power*100)}% Power')
        
        plt.xlabel('Sample Size (N)')
        plt.ylabel('Statistical Power')
        plt.title(f'Power Curves (α = {alpha})')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.ylim(0, 1)
        
        plt.tight_layout()
        plt.show()
    
    def bonferroni_correction(self, alpha: float, num_tests: int) -> float:
        """
        Calculate Bonferroni-corrected alpha level
        
        Parameters:
        -----------
        alpha : float
            Original alpha level
        num_tests : int
            Number of statistical tests
            
        Returns:
        --------
        float
            Bonferroni-corrected alpha level
        """
        return alpha / num_tests
    
    def analyze_study(self, target_r: float, n: Optional[int] = None, 
                     power: float = 0.8, alpha: float = 0.05, 
                     num_tests: int = 1, show_plot: bool = True) -> dict:
        """
        Comprehensive analysis for a study design
        
        Parameters:
        -----------
        target_r : float
            Target/expected correlation coefficient
        n : int, optional
            Sample size (if None, will calculate required N)
        power : float, default 0.8
            Desired statistical power
        alpha : float, default 0.05
            Significance level
        num_tests : int, default 1
            Number of statistical tests (for multiple comparisons)
        show_plot : bool, default True
            Whether to show power curve plot
            
        Returns:
        --------
        dict
            Comprehensive results dictionary
        """
        # Multiple comparisons correction
        alpha_corrected = self.bonferroni_correction(alpha, num_tests) if num_tests > 1 else alpha
        
        # Calculate sample size if not provided
        if n is None:
            n_required = self.sample_size(target_r, power, alpha)
            n_required_corrected = self.sample_size(target_r, power, alpha_corrected)
        else:
            n_required = n
            n_required_corrected = n
        
        # Calculate power
        power_uncorrected = self.power(n_required, target_r, alpha)
        power_corrected = self.power(n_required_corrected, target_r, alpha_corrected)
        
        # Minimum detectable effects
        min_effect_uncorrected = self.minimum_effect(n_required, power, alpha)
        min_effect_corrected = self.minimum_effect(n_required_corrected, power, alpha_corrected)
        
        results = {
            'target_correlation': target_r,
            'sample_size': n_required,
            'alpha_uncorrected': alpha,
            'alpha_corrected': alpha_corrected,
            'num_tests': num_tests,
            'power_uncorrected': power_uncorrected,
            'power_corrected': power_corrected,
            'min_detectable_uncorrected': min_effect_uncorrected,
            'min_detectable_corrected': min_effect_corrected,
            'sample_size_corrected': n_required_corrected
        }
        
        # Print summary
        print(f"CORRELATION POWER ANALYSIS RESULTS")
        print(f"=" * 40)
        print(f"Target correlation: r = {target_r:.3f}")
        print(f"Number of tests: {num_tests}")
        print(f"")
        print(f"SAMPLE SIZE REQUIREMENTS:")
        print(f"For 80% power (α = {alpha:.3f}): N = {n_required}")
        if num_tests > 1:
            print(f"For 80% power (α = {alpha_corrected:.4f}, Bonferroni): N = {n_required_corrected}")
        print(f"")
        print(f"POWER WITH N = {n_required}:")
        print(f"Uncorrected: {power_uncorrected:.3f}")
        if num_tests > 1:
            print(f"Bonferroni corrected: {power_corrected:.3f}")
        print(f"")
        print(f"MINIMUM DETECTABLE EFFECTS (80% power):")
        print(f"Uncorrected: r ≥ {min_effect_uncorrected:.3f}")
        if num_tests > 1:
            print(f"Bonferroni corrected: r ≥ {min_effect_corrected:.3f}")
        
        # Plot if requested
        if show_plot:
            effect_sizes = [target_r * 0.7, target_r, target_r * 1.3]
            self.plot_power_curves(effect_sizes, alpha=alpha)
        
        return results

# Example usage and quick functions for common use cases
def quick_sample_size(r: float, power: float = 0.8, alpha: float = 0.05) -> int:
    """Quick function to calculate required sample size"""
    analyzer = CorrelationPowerAnalysis()
    return analyzer.sample_size(r, power, alpha)

def quick_power(n: int, r: float, alpha: float = 0.05) -> float:
    """Quick function to calculate power"""
    analyzer = CorrelationPowerAnalysis()
    return analyzer.power(n, r, alpha)

def quick_min_effect(n: int, power: float = 0.8, alpha: float = 0.05) -> float:
    """Quick function to calculate minimum detectable effect"""
    analyzer = CorrelationPowerAnalysis()
    return analyzer.minimum_effect(n, power, alpha)

# Example usage
if __name__ == "__main__":
    # Create analyzer
    analyzer = CorrelationPowerAnalysis()
    
    # Example 1: Basic power analysis
    print("Example 1: Sample size for r = 0.3, 80% power")
    n_needed = analyzer.sample_size(r=0.3, power=0.8)
    print(f"Required sample size: {n_needed}")
    print()
    
    # Example 2: Power with given sample size
    print("Example 2: Power with N = 100, r = 0.25")
    power_val = analyzer.power(n=100, r=0.25)
    print(f"Statistical power: {power_val:.3f}")
    print()
    
    # Example 3: Minimum detectable effect
    print("Example 3: Minimum detectable effect with N = 80")
    min_r = analyzer.minimum_effect(n=80)
    print(f"Minimum detectable correlation: {min_r:.3f}")
    print()
    
    # Example 4: Comprehensive study analysis
    print("Example 4: Comprehensive analysis for a study")
    results = analyzer.analyze_study(target_r=0.3, num_tests=3, show_plot=False)