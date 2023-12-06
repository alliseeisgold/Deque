/** @author yaishenka
    @date 05.01.2023 */
#include <gtest/gtest.h>
#include <random>
#include <type_traits>
#include "deque.hpp"
#include "utils.hpp"
#include "memory_utils.hpp"

size_t MemoryManager::type_new_allocated = 0;
size_t MemoryManager::type_new_deleted = 0;
size_t MemoryManager::allocator_allocated = 0;
size_t MemoryManager::allocator_deallocated = 0;
size_t MemoryManager::allocator_constructed = 0;
size_t MemoryManager::allocator_destroyed = 0;

template <typename T, bool PropagateOnConstruct, bool PropagateOnAssign>
size_t
        WhimsicalAllocator<T, PropagateOnConstruct, PropagateOnAssign>::counter = 0;

size_t Accountant::ctor_calls = 0;
size_t Accountant::dtor_calls = 0;

bool ThrowingAccountant::need_throw = false;

void SetupTest() {
    MemoryManager::type_new_allocated = 0;
    MemoryManager::type_new_deleted = 0;
    MemoryManager::allocator_allocated = 0;
    MemoryManager::allocator_deallocated = 0;
    MemoryManager::allocator_constructed = 0;
    MemoryManager::allocator_destroyed = 0;
}

struct NotDefaultConstructible {
    NotDefaultConstructible() = delete;
    NotDefaultConstructible(int data) : data(data) {}
    int data;

    auto operator<=>(const NotDefaultConstructible&) const = default;
};

struct ThrowStruct {
    ThrowStruct(int value, bool throw_in_assignment, bool throw_in_copy) :
            value(value),
            throw_in_assignment(throw_in_assignment),
            throw_in_copy(throw_in_copy) {}

    ThrowStruct(const ThrowStruct& s) {
        value = s.value;
        throw_in_assignment = s.throw_in_assignment;
        throw_in_copy = s.throw_in_copy;

        if (throw_in_copy) {
            throw 1;
        }
    }

    ThrowStruct& operator=(const ThrowStruct& s) {
        if (throw_in_assignment) {
            throw 1;
        }

        value = s.value;
        throw_in_assignment = s.throw_in_assignment;
        throw_in_copy = s.throw_in_copy;

        return *this;
    }

    auto operator<=>(const ThrowStruct&) const = default;

    int value;
    bool throw_in_assignment;
    bool throw_in_copy;
};

template <class Stack1, class Stack2>
bool CompareStacks(const Stack1& s1, const Stack2& s2) {
    for (size_t i = 0; i < s1.size(); ++i) {
        if (s1[i] != s2[i]) {
            return false;
        }
    }

    return true;
}

template <typename Iterator, typename T>
void IteratorTest() {
    using traits = std::iterator_traits<Iterator>;

    {
        auto test = std::is_same_v<std::remove_cv_t<typename traits::value_type>,
                std::remove_cv_t<T>>;
        EXPECT_TRUE(test);
    }

    {
        auto test = std::is_same_v<typename traits::pointer, T*>;
        EXPECT_TRUE(test);
    }

    {
        auto test = std::is_same_v<typename traits::reference, T&>;
        EXPECT_TRUE(test);
    }

    {
        auto test = std::is_same_v<typename traits::iterator_category,
                std::random_access_iterator_tag>;
        EXPECT_TRUE(test);
    }

    {
        auto test = std::is_same_v<decltype(std::declval<Iterator>()++), Iterator>;
        EXPECT_TRUE(test);
    }

    {
        auto test = std::is_same_v<decltype(++std::declval<Iterator>()), Iterator&>;
        EXPECT_TRUE(test);
    }

    {
        auto
                test = std::is_same_v<decltype(std::declval<Iterator>() + 5), Iterator>;
        EXPECT_TRUE(test);
    }

    {
        auto test =
                std::is_same_v<decltype(std::declval<Iterator>() += 5), Iterator&>;
        EXPECT_TRUE(test);
    }

    {
        auto test = std::is_same_v<decltype(std::declval<Iterator>()
                                            - std::declval<Iterator>()), typename traits::difference_type>;
        EXPECT_TRUE(test);
    }

    {
        auto test = std::is_same_v<decltype(*std::declval<Iterator>()), T&>;
        EXPECT_TRUE(test);
    }

    {
        auto test = std::is_same_v<decltype(std::declval<Iterator>()
                                            < std::declval<Iterator>()), bool>;
        EXPECT_TRUE(test);
    }

    {
        auto test = std::is_same_v<decltype(std::declval<Iterator>()
                                            <= std::declval<Iterator>()), bool>;
        EXPECT_TRUE(test);
    }

    {
        auto test = std::is_same_v<decltype(std::declval<Iterator>()
                                            > std::declval<Iterator>()), bool>;
        EXPECT_TRUE(test);
    }

    {
        auto test = std::is_same_v<decltype(std::declval<Iterator>()
                                            >= std::declval<Iterator>()), bool>;
        EXPECT_TRUE(test);
    }

    {
        auto test = std::is_same_v<decltype(std::declval<Iterator>()
                                            == std::declval<Iterator>()), bool>;
        EXPECT_TRUE(test);
    }

    {
        auto test = std::is_same_v<decltype(std::declval<Iterator>()
                                            != std::declval<Iterator>()), bool>;
        EXPECT_TRUE(test);
    }
}

TEST(DequeConstructors, Default) {
    Deque<int> defaulted;
    EXPECT_EQ(defaulted.size(), 0);

    Deque<NotDefaultConstructible> without_default;
    EXPECT_EQ(without_default.size(), 0);
}

TEST(DequeConstructors, CopyEmpty) {
    Deque<NotDefaultConstructible> without_default;
    Deque<NotDefaultConstructible> copy = without_default;
    std::ignore = copy;
    EXPECT_EQ(copy.size(), 0);
}

TEST(DequeConstructors, WithSize) {
    size_t size = 17;
    int value = 14;

    {
        Deque<int> simple(size);
        EXPECT_EQ(simple.size(), size);
        EXPECT_TRUE(std::all_of(simple.begin(),
                                simple.end(),
                                [](int item) { return item == 0; }));
    }

    {
        Deque<NotDefaultConstructible> less_simple(size, value);
        EXPECT_EQ(less_simple.size(), size);
        EXPECT_TRUE(std::all_of(less_simple.begin(),
                                less_simple.end(),
                                [value = value](const auto& item) {
                                    return item.data == value;
                                }));
    }

}

TEST(DequeOperators, Assignment) {
    Deque<int> first(10, 10);
    Deque<int> second(9, 9);
    first = second;

    EXPECT_EQ(first.size(), 9);
    EXPECT_EQ(first.size(), second.size());
    EXPECT_TRUE(CompareStacks(first, second));
}

TEST(DequeConstructors, StaticAsserts) {
    using T1 = int;
    using T2 = NotDefaultConstructible;

    EXPECT_TRUE(std::is_default_constructible_v<Deque<T1>>);
    EXPECT_TRUE(std::is_default_constructible_v<Deque<T2>>);

    EXPECT_TRUE(std::is_copy_constructible_v<Deque<T1>>);
    EXPECT_TRUE(std::is_copy_constructible_v<Deque<T2>>);

    {
        auto test = std::is_constructible_v<Deque<T1>, size_t>;
        EXPECT_TRUE(test);
    }

    {
        auto test = std::is_constructible_v<Deque<T1>, size_t, const T1&>;
        EXPECT_TRUE(test);
    }

    {
        auto test = std::is_constructible_v<Deque<T2>, size_t, const T2&>;
        EXPECT_TRUE(test);
    }

    EXPECT_TRUE(std::is_copy_assignable_v<Deque<T1>>);
    EXPECT_TRUE(std::is_copy_assignable_v<Deque<T2>>);
}

TEST(DequeAccess, SquareBrackets) {
    Deque<size_t> defaulted(1300, 43);

    EXPECT_EQ(defaulted[0], defaulted[1280]);
    EXPECT_EQ(defaulted[0], 43);
}

TEST(DequeAccess, OperatorAt) {
    Deque<size_t> defaulted(1300, 43);

    EXPECT_EQ(defaulted.at(0), defaulted.at(1280));
    EXPECT_EQ(defaulted.at(0), 43);

    EXPECT_THROW(defaulted.at(size_t(-1)), std::out_of_range);
    EXPECT_THROW(defaulted.at(1300), std::out_of_range);
}

TEST(DequeAccess, StaticAsserts) {
    Deque<size_t> defaulted;
    const Deque<size_t> constant;

    {
        auto test = std::is_same_v<decltype(defaulted[0]), size_t&>;
        EXPECT_TRUE(test);
    }

    {
        auto test = std::is_same_v<decltype(defaulted.at(0)), size_t&>;
        EXPECT_TRUE(test);
    }

    {
        auto test = std::is_same_v<decltype(constant[0]), const size_t&>;
        EXPECT_TRUE(test);
    }

    {
        auto test = std::is_same_v<decltype(constant.at(0)), const size_t&>;
        EXPECT_TRUE(test);
    }

    EXPECT_FALSE(noexcept(defaulted.at(0)));
}

TEST(DequeItetators, StaticAsserts) {
    IteratorTest<Deque<int>::iterator, int>();
    IteratorTest<decltype(std::declval<Deque<int>>().rbegin()), int>();
    IteratorTest<decltype(std::declval<Deque<int>>().cbegin()), const int>();
}

TEST(DequtIterators, Arithmetic) {
    {
        Deque<int> empty;

        EXPECT_EQ(empty.end() - empty.begin(), 0);
        EXPECT_EQ(empty.begin() + 0, empty.end());
        EXPECT_EQ(empty.end() - 0, empty.begin());

        EXPECT_EQ(empty.rend() - empty.rbegin(), 0);
        EXPECT_EQ(empty.rbegin() + 0, empty.rend());
        EXPECT_EQ(empty.rend() - 0, empty.rbegin());

        EXPECT_EQ(empty.cend() - empty.cbegin(), 0);
        EXPECT_EQ(empty.cbegin() + 0, empty.cend());
        EXPECT_EQ(empty.cend() - 0, empty.cbegin());
    }

    {
        Deque<int> one(1);
        auto iter = one.end();

        EXPECT_EQ(--iter, one.begin());
        EXPECT_EQ(iter++, one.begin());
    }

    {
        Deque<int> d(1000, 3);

        EXPECT_EQ(size_t(d.end() - d.begin()), d.size());
        EXPECT_EQ(d.begin() + d.size(), d.end());
        EXPECT_EQ(d.end() - d.size(), d.begin());
    }
}

TEST(DequtIterators, Comparison) {
    Deque<int> d(1000, 3);

    EXPECT_TRUE(d.end() > d.begin());
    EXPECT_TRUE(d.cend() > d.cbegin());
    EXPECT_TRUE(d.rend() > d.rbegin());
}

TEST(DequeIterators, Algorithm) {
    Deque<int> d(1000, 3);

    std::iota(d.begin(), d.end(), 13);
    std::mt19937 g(31415);
    std::shuffle(d.begin(), d.end(), g);
    std::sort(d.rbegin(), d.rbegin() + 500);
    std::reverse(d.begin(), d.end());
    auto sorted_border = std::is_sorted_until(d.begin(), d.end());

    EXPECT_EQ(sorted_border - d.begin(), 500);
}

TEST(DequeModification, PushAndPop) {
    Deque<NotDefaultConstructible> d(10000, {1});
    int start_size = static_cast<int>(d.size());

    auto middle_iter = d.begin() + (start_size / 2); // 5000
    auto& middle_element = *middle_iter;
    auto begin = d.begin();
    auto end = d.rbegin();

    auto middle2_iter = middle_iter + 2000; // 7000

    // remove 400 elements
    for (size_t i = 0; i < 400; ++i) {
        d.pop_back();
    }

    // begin and middle iterators are still valid
    EXPECT_TRUE(begin->data == 1);
    EXPECT_TRUE(middle_iter->data == 1);
    EXPECT_TRUE(middle_element.data == 1);
    EXPECT_TRUE(middle2_iter->data == 1);

    end = d.rbegin();

    // 800 elemets removed in total
    for (size_t i = 0; i < 400; ++i) {
        d.pop_front();
    }

    // and and middle iterators are still valid
    EXPECT_TRUE(end->data == 1);
    EXPECT_TRUE(middle_iter->data == 1);
    EXPECT_TRUE(middle_element.data == 1);
    EXPECT_TRUE(middle2_iter->data == 1);

    // removed 9980 items in total
    for (size_t i = 0; i < 4590; ++i) {
        d.pop_front();
        d.pop_back();
    }

    EXPECT_TRUE(d.size() == 20);
    EXPECT_TRUE(middle_element.data == 1);
    EXPECT_TRUE(middle_iter->data == 1 && middle_iter->data == 1);
    EXPECT_TRUE(std::all_of(d.begin(),
                            d.end(),
                            [](const auto& item) { return item.data == 1; }));

    auto& begin_ref = *d.begin();
    auto& end_ref = *d.rbegin();

    for (size_t i = 0; i < 5500; ++i) {
        d.push_back({2});
        d.push_front({2});
    }

    EXPECT_TRUE((begin_ref).data == 1);
    EXPECT_TRUE((end_ref).data == 1);
    EXPECT_TRUE(d.begin()->data == 2);
    EXPECT_TRUE(d.size() == 5500 * 2 + 20);
}

TEST(DequeModification, InsertAndErase) {
    Deque<NotDefaultConstructible> d(10000, {1});
    auto start_size = d.size();

    d.insert(d.begin() + static_cast<int>(start_size) / 2,
             NotDefaultConstructible{2});
    EXPECT_TRUE(d.size() == start_size + 1);
    d.erase(d.begin() + static_cast<int>(start_size) / 2 - 1);
    EXPECT_TRUE(d.size() == start_size);

    EXPECT_TRUE(size_t(std::count(d.begin(), d.end(), NotDefaultConstructible{1}))
                == start_size - 1);
    EXPECT_TRUE(std::count(d.begin(), d.end(), NotDefaultConstructible{2}) == 1);

    Deque<NotDefaultConstructible> copy;
    for (const auto& item: d) {
        copy.insert(copy.end(), item);
    }

    EXPECT_TRUE(d.size() == copy.size());
    EXPECT_TRUE(std::equal(d.begin(), d.end(), copy.begin()));
}

TEST(Deque, Throw) {
    {
        EXPECT_THROW(Deque<ThrowStruct> d(10, ThrowStruct(0, false, true)), int);
    }

    {
        Deque<ThrowStruct> d(10, ThrowStruct(10, true, false));

        EXPECT_THROW(d[0] = ThrowStruct(1, false, false), int);
        EXPECT_EQ(d.size(), 10);
        EXPECT_EQ(d[0].value, 10);
    }

    {
        Deque<ThrowStruct> d(1, ThrowStruct(10, false, false));
        EXPECT_THROW(d.push_back(ThrowStruct(1, false, true)), int);
        EXPECT_EQ(d.size(), 1);
    }
}

/*
 * Pt2 tests
*/

TEST(DequeConstructors, DefaultConstructPt2) {
    SetupTest();
    {
        Deque<TypeWithFancyNewDeleteOperators> d;

        ASSERT_EQ(d.size(), 0);
        ASSERT_TRUE(d.empty());
        ASSERT_EQ(MemoryManager::type_new_allocated, 0);
    }
    ASSERT_EQ(MemoryManager::type_new_deleted, 0);
}

TEST(DequeConstructors, ConstructFromSizePt2) {
    SetupTest();
    {
        constexpr size_t kSize = 10;
        Deque<TypeWithFancyNewDeleteOperators> d(kSize);
        ASSERT_EQ(d.size(), kSize);
        ASSERT_EQ(MemoryManager::type_new_allocated, 0);
    }
    ASSERT_EQ(MemoryManager::type_new_deleted, 0);
}

TEST(DequeConstructors, ConstructFromSizeAndAlloc) {
    SetupTest();
    AllocatorWithCount<TypeWithFancyNewDeleteOperators> alloc;
    constexpr size_t kSize = 10;
    {
        Deque<TypeWithFancyNewDeleteOperators,
        AllocatorWithCount<TypeWithFancyNewDeleteOperators>> d(kSize, alloc);
        ASSERT_EQ(d.size(), kSize);
        ASSERT_EQ(MemoryManager::type_new_allocated, 0);
        ASSERT_EQ(MemoryManager::type_new_deleted, 0);
        ASSERT_TRUE(MemoryManager::allocator_allocated != 0);
        ASSERT_TRUE(MemoryManager::allocator_constructed >= kSize);
        ASSERT_EQ(alloc.allocator_allocated, 0);
        ASSERT_EQ(alloc.allocator_constructed, 0);
    }
    ASSERT_TRUE(MemoryManager::allocator_deallocated != 0);
    ASSERT_TRUE(MemoryManager::allocator_destroyed >= kSize);
    ASSERT_EQ(alloc.allocator_destroyed, 0);
    ASSERT_EQ(alloc.allocator_deallocated, 0);
}

TEST(Construct, ConstructFromSizeAndAllocAndDefaultValue) {
    SetupTest();
    constexpr size_t kSize = 10;
    AllocatorWithCount<TypeWithFancyNewDeleteOperators> alloc;
    TypeWithFancyNewDeleteOperators default_value(1);
    Deque<TypeWithFancyNewDeleteOperators,
    AllocatorWithCount<TypeWithFancyNewDeleteOperators>>
            d(kSize, default_value, alloc);
    ASSERT_EQ(d.size(), kSize);
    ASSERT_EQ(MemoryManager::type_new_allocated, 0);
    ASSERT_EQ(MemoryManager::type_new_deleted, 0);
    ASSERT_TRUE(MemoryManager::allocator_allocated != 0);
    ASSERT_TRUE(MemoryManager::allocator_constructed >= kSize);
    ASSERT_EQ(alloc.allocator_allocated, 0);
    ASSERT_EQ(alloc.allocator_constructed, 0);

    for (auto it = d.begin(); it != d.end(); ++it) {
        ASSERT_EQ(it->value, default_value.value);
    }
}

TEST(Construct, Copy) {
    SetupTest();
    const size_t size = 5;
    Deque<TypeWithCounts, AllocatorWithCount<TypeWithCounts>> d1 = {1, 2, 3, 4, 5};
    Deque<TypeWithCounts, AllocatorWithCount<TypeWithCounts>> d2(d1);

    ASSERT_EQ(d2.size(), size);
    ASSERT_TRUE(MemoryManager::allocator_allocated != 0);
    ASSERT_TRUE(MemoryManager::allocator_constructed >= size * 2);
    ASSERT_TRUE(CompareStacks(d1, d2));
    for (auto& value: d1) {
        ASSERT_EQ(*value.copy_c, 2);
    }
}

TEST(Propagate, PropagateOnConstructAndPropagateOnAssign) {
    SetupTest();
    Deque<int, WhimsicalAllocator<int, true, true>> d;

    d.push_back(1);
    d.push_back(2);

    auto copy = d;
    assert(copy.get_allocator() != d.get_allocator());

    d = copy;
    assert(copy.get_allocator() == d.get_allocator());
}

TEST(Propagate, NotPropagateOnConstructAndNotPropagateOnAssign) {
    SetupTest();
    Deque<int, WhimsicalAllocator<int, false, false>> d;

    d.push_back(1);
    d.push_back(2);

    auto copy = d;
    assert(copy.get_allocator() == d.get_allocator());

    d = copy;
    assert(copy.get_allocator() == d.get_allocator());
}

TEST(Propagate, PropagateOnConstructAndNotPropagateOnAssign) {
    SetupTest();
    Deque<int, WhimsicalAllocator<int, true, false>> d;

    d.push_back(1);
    d.push_back(2);

    auto copy = d;
    assert(copy.get_allocator() != d.get_allocator());

    d = copy;
    assert(copy.get_allocator() != d.get_allocator());
}

TEST(Deque, TestAccountant) {
    Accountant::reset();
    {
        Deque<Accountant> d(5);
        ASSERT_EQ(d.size(), 5);
        ASSERT_EQ(Accountant::ctor_calls, 5);

        Deque<Accountant> another = d;
        ASSERT_EQ(another.size(), 5);
        ASSERT_EQ(Accountant::ctor_calls, 10);
        ASSERT_EQ(Accountant::dtor_calls, 0);

        another.pop_back();
        another.pop_front();
        ASSERT_EQ(Accountant::dtor_calls, 2);

        d = another; // dtor_calls += 5, ctor_calls += 3
        ASSERT_EQ(another.size(), 3);
        ASSERT_EQ(d.size(), 3);

        ASSERT_EQ(Accountant::ctor_calls, 13);
        ASSERT_EQ(Accountant::dtor_calls, 7);
    } // dtor_calls += 6

    ASSERT_EQ(Accountant::ctor_calls, 13);
    ASSERT_EQ(Accountant::dtor_calls, 13);
}

TEST(Deque, ExceptionSafety) {
    Accountant::reset();

    ThrowingAccountant::need_throw = true;

    try {
        Deque<ThrowingAccountant> d(8);
    } catch (...) {
        ASSERT_EQ(Accountant::ctor_calls, 4);
        ASSERT_EQ(Accountant::dtor_calls, 4);
    }

    ThrowingAccountant::need_throw = false;
    Deque<ThrowingAccountant> d(8);

    Deque<ThrowingAccountant> d2;
    for (int i = 0; i < 13; ++i) {
        d2.push_back(i);
    }

    Accountant::reset();
    ThrowingAccountant::need_throw = true;

    try {
        auto d3 = d2;
    } catch (...) {
        ASSERT_EQ(Accountant::ctor_calls, 4);
        ASSERT_EQ(Accountant::dtor_calls, 4);
    }

    Accountant::reset();

    try {
        d = d2;
    } catch (...) {
        ASSERT_EQ(Accountant::ctor_calls, 4);
        ASSERT_EQ(Accountant::dtor_calls, 4);

        // Actually it may not be 8 (although de facto it is), but the only thing we can demand here
        // is the abscence of memory leaks
        //
        //assert(lst.size() == 8);
    }
}

TEST(DequeConstructors, OnlyMovable) {
    SetupTest();
    Deque<OnlyMovable> d;
    // NO LINT
    ASSERT_EQ(d.size(), 0);
    d.emplace_back(0);
    ASSERT_EQ(d.size(), 1);
    OnlyMovable om(0);
    d.push_back(std::move(om));
    ASSERT_EQ(d.size(), 2);
}

TEST(DequeConstructors, FromRvalueDeque) {
    SetupTest();
    const size_t size = 5;
    Deque<TypeWithCounts, AllocatorWithCount<TypeWithCounts>> d1 = {1, 2, 3, 4, 5};
    Deque<TypeWithCounts, AllocatorWithCount<TypeWithCounts>> d2(std::move(d1));

    ASSERT_EQ(d1.size(), 0);
    ASSERT_EQ(d2.size(), size);
    ASSERT_TRUE(MemoryManager::allocator_allocated != 0);
    ASSERT_TRUE(MemoryManager::allocator_constructed >= size);
    for (auto& value : d1) {
        ASSERT_EQ(*value.copy_c, 1);
        ASSERT_EQ(*value.move_c, 1);
    }

}

TEST(Deque, EmplaceBackPushBackRValue) {
    TypeWithCounts t(0);
    Deque<TypeWithCounts> d;
    d.emplace_back(0);
    ASSERT_EQ(*d.begin()->int_c, 1);
    ASSERT_EQ(*d.begin()->move_c, 0);
    ASSERT_EQ(*d.begin()->copy_c, 0);

    d.push_front(std::move(t));
    ASSERT_EQ(*d.begin()->move_c, 1);
    ASSERT_EQ(*d.begin()->copy_c, 0);
}


int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}