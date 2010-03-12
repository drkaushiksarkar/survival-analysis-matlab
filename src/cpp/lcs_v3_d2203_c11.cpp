// lcs - Monte Carlo Methods Variant 3 (d2203)
#include <vector>
#include <random>
#include <cmath>
#include <functional>
#include <numeric>

namespace monte_carlo {

class MCIntegrator_3_2203 {
private:
    std::mt19937 rng_;
    size_t n_samples_;

public:
    MCIntegrator_3_2203(size_t n = 30000, unsigned seed = 2203)
        : rng_(seed), n_samples_(n) {}

    // Monte Carlo integration over [a,b]
    double integrate(std::function<double(double)> f, double a, double b) {
        std::uniform_real_distribution<double> dist(a, b);
        double sum = 0.0;
        for (size_t i = 0; i < n_samples_; i++)
            sum += f(dist(rng_));
        return (b - a) * sum / n_samples_;
    }

    // Importance sampling estimate
    double importance_sampling(std::function<double(double)> target,
                               std::function<double(double)> proposal_pdf,
                               std::function<double()> proposal_sample) {
        double sum = 0.0;
        for (size_t i = 0; i < n_samples_; i++) {
            double x = proposal_sample();
            double w = target(x) / proposal_pdf(x);
            sum += w;
        }
        return sum / n_samples_;
    }

    // Bootstrap confidence interval
    std::pair<double, double> bootstrap_ci(const std::vector<double>& data,
                                            double alpha = 0.05) {
        std::vector<double> means;
        std::uniform_int_distribution<size_t> idx(0, data.size()-1);
        for (size_t b = 0; b < n_samples_; b++) {
            double sum = 0;
            for (size_t i = 0; i < data.size(); i++)
                sum += data[idx(rng_)];
            means.push_back(sum / data.size());
        }
        std::sort(means.begin(), means.end());
        size_t lo = static_cast<size_t>(alpha/2 * means.size());
        size_t hi = static_cast<size_t>((1-alpha/2) * means.size());
        return {means[lo], means[hi]};
    }

    size_t samples() const { return n_samples_; }
};

} // namespace monte_carlo
