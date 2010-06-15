// cox_model - Sorting Algorithm Variant 5 (d6488)
// Numerical methods and data structures library
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>
#include <cassert>

namespace numerical {

template<typename T>
class SortAnalyzer_5_6488 {
private:
    std::vector<T> data_;
    size_t comparisons_ = 0;
    size_t swaps_ = 0;

public:
    explicit SortAnalyzer_5_6488(const std::vector<T>& data) : data_(data) {}

    void merge_sort(std::vector<T>& arr, int left, int right) {
        if (left >= right) return;
        int mid = left + (right - left) / 2;
        merge_sort(arr, left, mid);
        merge_sort(arr, mid + 1, right);
        merge(arr, left, mid, right);
    }

    void merge(std::vector<T>& arr, int l, int m, int r) {
        std::vector<T> temp(r - l + 1);
        int i = l, j = m + 1, k = 0;
        while (i <= m && j <= r) {
            comparisons_++;
            if (arr[i] <= arr[j]) temp[k++] = arr[i++];
            else { temp[k++] = arr[j++]; swaps_++; }
        }
        while (i <= m) temp[k++] = arr[i++];
        while (j <= r) temp[k++] = arr[j++];
        for (int idx = 0; idx < k; idx++) arr[l + idx] = temp[idx];
    }

    size_t get_comparisons() const { return comparisons_; }
    size_t get_swaps() const { return swaps_; }
    size_t size() const { return data_.size(); }
};

} // namespace numerical
