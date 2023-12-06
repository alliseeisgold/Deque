#pragma once

#include <iostream>
#include <vector>

const int kBucketSize = 25;
const int kPtrsArraySize = 10;

struct Point {
    int pointer;
    int index;

    Point(int new_pointer, int new_index)
            : pointer(new_pointer), index(new_index) {}
};

template<typename T, typename Allocator = std::allocator<T>>

class Deque {
public:
    // Iterator
    template<bool IsConst>
    class Iterator {
    public:
        using value_type = typename std::conditional<IsConst, const T, T>::type;
        using pointer = typename std::conditional<IsConst, const T *, T *>::type;
        using reference = typename std::conditional<IsConst, const T &, T &>::type;
        using iterator_category = std::random_access_iterator_tag;
        using difference_type = std::ptrdiff_t;

        Iterator() = default;

        Iterator(Point current) : current_(current) {}

        Iterator(Point current, Deque<T, Allocator> *eternal)
                : current_(current), eternal_(eternal) {}

        void set_eternal(Deque<T, Allocator> *new_eternal) {
            eternal_ = new_eternal;
        }

        explicit operator Iterator<true>() const {
            return Iterator<true>(current_);
        }

        void set_pointer(int new_pointer) { current_.pointer = new_pointer; }

        void set_index(int new_index) { current_.index = new_index; }

        size_t get_pointer() const { return static_cast<size_t>(current_.pointer); }

        size_t get_index() { return static_cast<size_t>(current_.index); }

        Iterator &operator++() {
            if (current_.index + 1 < kBucketSize) {
                current_.index++;
                return *this;
            }
            current_.pointer++;
            current_.index = 0;
            return *this;
        }

        Iterator operator++(int) {
            Iterator result = *this;
            this->operator++();
            return result;
        }

        Iterator &operator--() {
            if (current_.index > 0) {
                current_.index--;
                return *this;
            }
            // если current_.pointer == 0?
            current_.pointer--;
            current_.index = kBucketSize - 1;
            return *this;
        }

        Iterator operator--(int) {
            Iterator result = *this;
            this->operator--();
            return result;
        }

        Iterator &operator+=(size_t n) {
            int step = static_cast<int>(n);
            current_.pointer += (current_.index + step) / kBucketSize;
            current_.index = (current_.index + step) % kBucketSize;
            return *this;
        }

        Iterator &operator-=(size_t n) {
            int step = static_cast<int>(n);
            if (current_.index >= step) {
                current_.index -= step;
                return *this;
            }
            // need to think
            step -= (current_.index + 1);
            current_.pointer -= step / kBucketSize;
            current_.pointer -= static_cast<int>(step % kBucketSize > 0);
            current_.index = kBucketSize - step % kBucketSize - 1;
            return *this;
        }

        Iterator operator+(size_t n) const {
            Iterator copy = *this;
            copy += n;
            return copy;
        }

        Iterator operator-(size_t n) const {
            Iterator copy = *this;
            copy -= n;
            return copy;
        }

        bool operator<(const Iterator &other) const {
            return current_.pointer < other.current_.pointer ||
                   (current_.pointer == other.current_.pointer &&
                    current_.index < other.current_.index);
        }

        bool operator>(const Iterator &other) const { return other < *this; }

        bool operator==(const Iterator &other) const {
            return current_.pointer == other.current_.pointer &&
                   current_.index == other.current_.index;
        }

        bool operator!=(const Iterator &other) const { return !(*this == other); }

        bool operator<=(const Iterator &other) const {
            return *this < other || *this == other;
        }

        bool operator>=(const Iterator &other) const {
            return *this > other || *this == other;
        }

        difference_type operator-(const Iterator &other) const {
            difference_type result;
            if (current_.pointer == other.current_.pointer) {
                result = current_.index - other.current_.index;
            } else {
                result = current_.pointer * kBucketSize + current_.index;
                result -= other.current_.pointer * kBucketSize + other.current_.index;
            }
            return result;
        }

        reference operator*() {
            return *((eternal_->data_)[current_.pointer] + current_.index);
        }

        pointer operator->() const {
            return ((eternal_->data_)[current_.pointer] + current_.index);
        }

    private:
        Point current_;
        Deque<T, Allocator> *eternal_;
    };

    using iterator = Iterator<false>;
    using const_iterator = Iterator<true>;
    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    iterator begin() { return start_; }

    iterator end() { return end_; }

    const_iterator begin() const { return start_; }

    const_iterator end() const { return end_; }

    const_iterator cbegin() const { return Iterator<true>(start_); }

    const_iterator cend() const { return Iterator<true>(end_); }

    reverse_iterator rbegin() { return std::make_reverse_iterator(end()); }

    reverse_iterator rend() { return std::make_reverse_iterator(begin()); }

    const_reverse_iterator rbegin() const {
        return std::make_reverse_iterator(end());
    }

    const_reverse_iterator rend() const {
        return std::make_reverse_iterator(begin());
    }

    using alloc_traits = std::allocator_traits<Allocator>;

    // Constructors
    Deque(const Allocator &allocator = Allocator())
            : alloc_(allocator),
              size_(kPtrsArraySize),
              start_(Point(0, (kBucketSize - 1) / 2), this),
              end_(Point(0, (kBucketSize - 1) / 2), this),
              counter_(0),
              data_(kPtrsArraySize) {
        allocate_it();
    }

    Deque(const Deque &other)
            : alloc_(
            alloc_traits::select_on_container_copy_construction(other.alloc_)),
              counter_(other.counter_),
              size_(other.size_),
              start_(Point(0, (kBucketSize - 1) / 2), this),
              end_(Point(0, (kBucketSize - 1) / 2), this),
              data_(size_) {
        allocate_it();
        try {
            for (size_t idx = 0; idx < counter_; idx++) {
                alloc_traits::construct(
                        alloc_, data_[end_.get_pointer()] + end_.get_index(), other[idx]);
                end_++;
            }
        } catch (...) {
            remove_all();
            throw;
        }
    }

    Deque(size_t count, const Allocator &allocator = Allocator())
            : alloc_(allocator),
              counter_(count),
              size_(2 * count / kBucketSize + 2),
              start_(Point(0, (kBucketSize - 1) / 2), this),
              end_(Point(0, (kBucketSize - 1) / 2), this),
              data_(size_) {
        allocate_it();
        try {
            for (size_t idx = 0; idx < counter_; idx++) {
                alloc_traits::construct(alloc_,
                                        data_[end_.get_pointer()] + end_.get_index());
                end_++;
            }
        } catch (...) {
            remove_all();
            throw;
        }
    }

    Deque(size_t count, const T &value, const Allocator &allocator = Allocator())
            : alloc_(allocator),
              counter_(count),
              size_(2 * count / kBucketSize + 2),
              start_(Point(0, (kBucketSize - 1) / 2), this),
              end_(Point(0, (kBucketSize - 1) / 2), this),
              data_(std::vector<T *>(size_)) {
        allocate_it();
        try {
            for (size_t idx = 0; idx < counter_; idx++) {
                alloc_traits::construct(
                        alloc_, data_[end_.get_pointer()] + end_.get_index(), value);
                end_++;
            }
        } catch (...) {
            remove_all();
            throw;
        }
    }

    Deque(Deque &&other)
            : alloc_(other.alloc_),
              counter_(other.size()),
              size_(other.size_),
              start_(other.start_),
              end_(other.end_),
              data_(std::move(other.data_)) {
        start_.set_eternal(this);
        end_.set_eternal(this);
        other.start_ = other.end_ = {Point(0, (kBucketSize - 1) / 2), &other};
        other.counter_ = 0;
        other.size_ = kPtrsArraySize;
    }

    Deque(std::initializer_list<T> init, const Allocator &allocator = Allocator())
            : alloc_(allocator),
              counter_(init.size()),
              size_(2 * init.size() / kBucketSize + 2),
              start_(Point(0, (kBucketSize - 1) / 2), this),
              end_(Point(0, (kBucketSize - 1) / 2), this),
              data_(size_) {
        allocate_it();
        try {
            for (const auto &elems: init) {
                alloc_traits::construct(
                        alloc_, data_[end_.get_pointer()] + end_.get_index(), elems);
                end_++;
            }
        } catch (...) {
            remove_all();
            throw;
        }
    }

    // Destructor
    ~Deque() { remove_all(); }

    // Assignment operators
    Deque &operator=(const Deque &other) {
        Deque current(other);
        if (alloc_traits::propagate_on_container_copy_assignment::value &&
            current.alloc_ != other.alloc_) {
            current.alloc_ = other.alloc_;
        }
        std::swap(alloc_, current.alloc_);
        std::swap(size_, current.size_);
        std::swap(start_, current.start_);
        std::swap(end_, current.end_);
        start_.set_eternal(this);
        end_.set_eternal(this);
        std::swap(counter_, current.counter_);
        std::swap(data_, current.data_);
        return *this;
    }

    // Size and empty-check functions
    size_t size() const { return counter_; }

    bool empty() const { return counter_ == 0; }

    // Access operators
    T &operator[](size_t index) {
        Iterator<false> copy = start_;
        copy += index;
        return *(data_[copy.get_pointer()] + copy.get_index());
    }

    const T &operator[](size_t index) const {
        Iterator<false> copy = start_;
        copy += index;
        return *(data_[copy.get_pointer()] + copy.get_index());
    }

    T &at(size_t index) {
        if (index >= counter_) {
            throw std::out_of_range("Index out of range");
        }
        return Deque<T>::operator[](index);
    }

    const T &at(size_t index) const {
        if (index >= counter_) {
            throw std::out_of_range("Index out of range");
        }
        return Deque<T>::operator[](index);
    }

    Allocator get_allocator() { return alloc_; }

    // Modifiers
    // back
    void push_back(const T &value) {
        if (end_.get_pointer() == size_ && end_.get_index() == 0) {
            for (size_t _ = 0; _ < size_; _++) {
                data_.push_back(alloc_traits::allocate(alloc_, kBucketSize));
            }
            size_ += size_;
        }
        alloc_traits::construct(
                alloc_, data_[end_.get_pointer()] + end_.get_index(), value);
        end_++;
        counter_++;
    }

    void push_back(T &&value) {
        if (end_.get_pointer() == size_ && end_.get_index() == 0) {
            for (size_t _ = 0; _ < size_; _++) {
                data_.push_back(alloc_traits::allocate(alloc_, kBucketSize));
            }
            size_ += size_;
        }
        alloc_traits::construct(alloc_,
                                data_[end_.get_pointer()] + end_.get_index(),
                                std::forward<T>(value));
        end_++;
        counter_++;
    }

    template<class... Args>
    void emplace_back(Args &&... values) {
        if (end_.get_pointer() == size_ && end_.get_index() == 0) {
            for (size_t _ = 0; _ < size_; _++) {
                data_.push_back(alloc_traits::allocate(alloc_, kBucketSize));
            }
            size_ += size_;
        }
        alloc_traits::construct(alloc_,
                                data_[end_.get_pointer()] + end_.get_index(),
                                std::forward<Args>(values)...);
        end_++;
        counter_++;
    }

    void pop_back() {
        alloc_traits::destroy(alloc_, data_[end_.get_pointer()] + end_.get_index());
        end_--;
        counter_--;
    }

    // front
    void push_front(const T &value) {
        if (start_.get_pointer() == 0 && start_.get_index() == 0) {
            for (size_t _ = 0; _ < size_; _++) {
                data_.push_back(alloc_traits::allocate(alloc_, kBucketSize));
            }
            start_.set_pointer(static_cast<int>(start_.get_pointer() + size_));
            end_.set_pointer(static_cast<int>(end_.get_pointer() + size_));
            size_ += size_;
        }
        start_--;
        alloc_traits::construct(
                alloc_, data_[start_.get_pointer()] + start_.get_index(), value);
        counter_++;
    }

    void push_front(T &&value) {
        if (start_.get_pointer() == 0 && start_.get_index() == 0) {
            for (size_t _ = 0; _ < size_; _++) {
                data_.push_back(alloc_traits::allocate(alloc_, kBucketSize));
            }
            start_.set_pointer(static_cast<int>(start_.get_pointer() + size_));
            end_.set_pointer(static_cast<int>(end_.get_pointer() + size_));
            size_ += size_;
        }
        start_--;
        alloc_traits::construct(alloc_,
                                data_[start_.get_pointer()] + start_.get_index(),
                                std::forward<T>(value));
        counter_++;
    }

    template<class... Args>
    void emplace_front(Args &&... values) {
        if (start_.get_pointer() == 0 && start_.get_index() == 0) {
            for (size_t _ = 0; _ < size_; _++) {
                data_.push_back(alloc_traits::allocate(alloc_, kBucketSize));
            }
            start_.set_pointer(static_cast<int>(start_.get_pointer() + size_));
            end_.set_pointer(static_cast<int>(end_.get_pointer() + size_));
            size_ += size_;
        }
        start_--;
        alloc_traits::construct(alloc_,
                                data_[start_.get_pointer()] + start_.get_index(),
                                std::forward<Args>(values)...);
        counter_++;
    }

    void pop_front() {
        alloc_traits::destroy(alloc_,
                              data_[start_.get_pointer()] + start_.get_index());
        start_++;
        counter_--;
    }

    // insert and erase
    void insert(iterator iterator1, const T &value) {
        if (iterator1 == begin()) {
            push_front(value);
        } else if (iterator1 == end()) {
            push_back(value);
        } else {
            auto iter = end() - 1;
            push_back(*iter);
            iter--;
            for (; iter > iterator1; --iter) {
                *iter = *(iter - 1);
            }
            *iterator1 = value;
        }
    }

    void erase(iterator iterator1) {
        if (iterator1 == begin()) {
            pop_front();
        } else if (iterator1 + 1 == end()) {
            pop_back();
        } else {
            for (auto iter = iterator1 + 1; iter < end(); ++iter) {
                *(iter - 1) = *iter;
            }
            pop_back();
        }
    }

private:
    // to avoid duplicates
    void allocate_it() {
        for (auto &ptr: data_) {
            ptr = alloc_traits::allocate(alloc_, kBucketSize);
        }
    }

    void remove_all() {
        for (auto &item: *this) {
            alloc_traits::destroy(alloc_, &item);
        }
        for (auto &ptr: data_) {
            alloc_traits::deallocate(alloc_, ptr, kBucketSize);
        }
        counter_ = 0;
        size_ = kPtrsArraySize;
    }

    Allocator alloc_;
    size_t size_;
    Iterator<false> start_;
    Iterator<false> end_;
    size_t counter_;
    std::vector<T *> data_;
};