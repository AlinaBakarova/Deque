#include <iostream>
#include <vector>
#include <algorithm>

using std::vector;

template <typename T, typename B>
class common_iterator;

template <typename T>
class Deque {
private:
    static const size_t max_size = 64;

    template <typename TD, typename B>
    friend class common_iterator;
    vector<T*> deque_;

    size_t first_arr = 0;
    size_t last_arr = 0;
    size_t front_ = 0;
    size_t back_ = 0;
    size_t size_ = 0;

public:
    Deque() {
        deque_.resize(1);
        deque_[0] = reinterpret_cast<T*>(new int8_t[sizeof(T) * max_size]);
    }

    Deque(const Deque &d)
            : first_arr(d.first_arr)
            , last_arr(d.last_arr)
            , front_(d.front_)
            , back_(d.back_)
            , size_ (d.size_)
    {
        deque_.resize(d.deque_.size());
        for (size_t i = first_arr; i <= last_arr; ++i) {
            deque_[i] = reinterpret_cast<T*>(new int8_t[max_size * sizeof(T)]);

            for (size_t j = 0; j < max_size; ++j) {
                if (i == first_arr && j < front_) continue;
                if (i == last_arr && j > back_) break;

                new(deque_[i] + j) T(d.deque_[i][j]);
            }
        }
    }

    Deque(size_t n, const T& value = T())
            : first_arr(0)
            , last_arr(std::max((n + max_size - 1) / max_size, size_t(1)) - 1)
            , front_(0)
            , back_((n - 1) % max_size)
            , size_(n)
    {
        deque_.resize(last_arr + 1);

        for (size_t y = first_arr; y <= last_arr; ++y) {
            deque_[y] = reinterpret_cast<T*>(new int8_t[max_size * sizeof(T)]);
            for (size_t x = 0; x < max_size; ++x) {
                if (y == last_arr && x > back_) break;
                new(deque_[y] + x) T(value);
            }
        }
    }

    Deque& operator=(const Deque& d) {
        Deque copy = Deque(d);
        std::swap(size_, copy.size_);
        std::swap(deque_, copy.deque_);
        std::swap(first_arr, copy.first_arr);
        std::swap(last_arr, copy.last_arr);
        std::swap(back_, copy.back_);
        std::swap(front_, copy.front_);
        return *this;
    }

    size_t size() const { return size_; }

    T& operator[](size_t i) {
        size_t position = first_arr * max_size + front_ + i;
        return deque_[position / max_size][position % max_size];
    }

    const T& operator[](size_t i) const {
        size_t position = first_arr * max_size + front_ + i;
        return deque_[position / max_size][position % max_size];
    }

    T& at(size_t i) {
        if (i >= size_) throw std::out_of_range("Error: out of range");
        return operator[](i);
    }

    const T& at(size_t i) const {
        if (i >= size_) throw std::out_of_range("Error: out of range");
        return operator[](i);
    }

    void push_back(const T& value) {
        if (size_ == 0) {
            new(deque_[last_arr] + back_) T(value);
            ++size_;
            return;
        }

        if (last_arr == deque_.size() - 1 && back_ == max_size - 1) {
            size_t prev_arr_size = deque_.size();
            deque_.resize(2 * deque_.size());
            for (size_t i = prev_arr_size; i < deque_.size(); ++i) {
                deque_[i] = reinterpret_cast<T*>(new int8_t[max_size * sizeof(T)]);
            }
        }

        ++back_;
        if (back_ == max_size) {
            ++last_arr;
            back_ = 0;
            if (!deque_[last_arr]) deque_[last_arr] = reinterpret_cast<T*>(new int8_t[max_size * sizeof(T)]);
        }

        new(deque_[last_arr] + back_) T(value);
        ++size_;
    }

    void pop_back() {
        (deque_[last_arr] + back_)->~T();
        if (back_ == 0) {
            delete deque_[last_arr];
            --last_arr;
            back_ = max_size;
        }
        --back_;
        --size_;
    }

    void push_front(const T &value) {
        if (size_ == 0) {
            new(deque_[first_arr] + front_) T(value);
            ++size_;
            return;
        }

        if (first_arr == 0 && front_ == 0) {
            size_t prev_size = deque_.size();
            vector<T*> new_deque(prev_size * 2);

            for (size_t i = 0; i < prev_size; i++) {
                new_deque[i + prev_size] = deque_[i];
            }

            first_arr = prev_size - 1;
            last_arr += prev_size;
            front_ = max_size - 1;

            deque_ = new_deque;
            deque_[first_arr] = reinterpret_cast<T*>(new int8_t[max_size * sizeof(T)]);
        } else if (first_arr != 0 && front_ == 0) {
            --first_arr;
            front_ = max_size - 1;

            deque_[first_arr] = reinterpret_cast<T*>(new int8_t[max_size * sizeof(T)]);
        } else {
            --front_;
        }

        new(deque_[first_arr] + front_) T(value);
        ++size_;
    }

    void pop_front() {
        (deque_[first_arr] + front_)->~T();
        --size_;
        ++front_;

        if (front_ == max_size) {
            delete deque_[first_arr];
            ++first_arr;
            front_ = 0;
        }
    }

    using iterator = common_iterator<T, Deque<T>*>;
    using reverse_iterator = std::reverse_iterator<common_iterator<T, Deque<T>*>>;
    using const_iterator = common_iterator<const T, const Deque<T>*>;
    using const_reverse_iterator = std::reverse_iterator<common_iterator<const T, const Deque<T>*>>;

    iterator begin() {
        return iterator(first_arr, front_, this);
    }

    iterator end() {
        auto copy = iterator(last_arr, back_, this);
        ++copy;
        return copy;
    }

    const_iterator begin() const{
        return const_iterator(first_arr, front_, this);
    }

    const_iterator end() const{
        auto copy = const_iterator(last_arr, back_, this);
        return ++copy;
    }

    const_iterator cbegin() const {
        return const_iterator(first_arr, front_, this);
    }

    const_iterator cend() const {
        auto copy = const_iterator(last_arr, back_, this);
        return ++copy;
    }

    reverse_iterator rbegin() { return std::reverse_iterator(end());}

    reverse_iterator rend() { return std::reverse_iterator(begin());}

    const_reverse_iterator rbegin() const { return std::reverse_iterator(cend());}

    const_reverse_iterator rend() const { return std::reverse_iterator(cbegin());}

    const_reverse_iterator crbegin() const { return std::reverse_iterator(cend());}

    const_reverse_iterator crend() const { return std::reverse_iterator(cbegin());}

    void insert(const iterator& it, const T& value) {
        push_back(value);
        for (auto i = --end(); i != it; i--) {
            std::iter_swap(i, i - 1);
        }
    }

    void erase(const iterator& it) {
        for (auto i = it; i != --end(); i++) {
            std::iter_swap(i, i + 1);
        }
        pop_back();
    }
};

template <typename T, typename deque_type>
class common_iterator : public std::iterator<std::random_access_iterator_tag, T> {
private:
public:
    size_t y_;
    size_t x_;
    deque_type deq_;

    static const size_t max_size = std::remove_pointer_t<deque_type>::max_size;

    common_iterator(size_t y, size_t x, const deque_type& deq): y_(y), x_(x), deq_(deq) {}

    common_iterator() = default;

    common_iterator& operator = (const common_iterator<T, deque_type>& it) {
        if (this == &it)
            return *this;
        x_ = it.x_;
        y_ = it.y_;
        deq_ = it.deq_;
        return *this;
    }

    common_iterator& operator--() {
        *this = *this - 1;
        return *this;
    }

    common_iterator operator--(int) {
        common_iterator res = *this;
        *this = *this - 1;
        return res;
    }

    common_iterator& operator++() {
        *this = *this + 1;
        return *this;
    }

    common_iterator operator++(int) {
        common_iterator res = *this;
        *this = *this + 1;
        return res;
    }

    common_iterator& operator-=(size_t p) {
        *this = *this + (-p);
        return *this;
    }
    common_iterator& operator+=(size_t p) {
        *this = *this + p;
        return *this;
    }

    T& operator*() const {
        return deq_->deque_[y_][x_];
    }

    size_t operator-(const common_iterator& it) const {
        return (y_ - it.y_ - 1) * max_size + (x_ + 1) + (max_size - 1 - it.x_);
    }

    common_iterator operator-(int n) const {
        int minus_block = (x_ - n) / max_size;
        int new_elem_index = (x_ - n) % max_size;
        if (new_elem_index < 0) {
            minus_block -= 1;
            new_elem_index += max_size;
        }
        return common_iterator(y_ + minus_block, new_elem_index, deq_);
    }

    T* operator->() {
        return &deq_->deque_[y_][x_];
    }

    bool operator==(const common_iterator& rit) const {
        return (rit.x_ == x_ && rit.y_ == y_);
    }

    bool operator<(const common_iterator& rit) const  {
        return (y_ == rit.y_ && x_ < rit.x_) || (y_ < rit.y_);
    }

    bool operator>(const common_iterator& rit) const {
        return !(*this <= rit);
    }

    bool operator>=(const common_iterator& rit) const {
        return !(*this < rit);
    }

    bool operator<=(const common_iterator& rit) const {
        return *this < rit || *this == rit;
    }

    bool operator!=(const common_iterator& rit) const {
        return !(*this == rit);
    }

    common_iterator operator+(size_t x) const {
        return common_iterator(y_ + (x_ + x) / max_size, (x_ + x) % max_size, deq_);
    }
};

template <typename T, typename deque_t>
common_iterator<T, deque_t> operator+(size_t p, const common_iterator<T, deque_t>& it) {
    return it + p;
}
