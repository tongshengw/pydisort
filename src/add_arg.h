#pragma once

// C/C++
#include <utility>

#ifdef __CUDACC__
#define DISPATCH_MACRO __host__ __device__
#else
#define DISPATCH_MACRO
#endif

#define ADD_ARG(T, name)                                                       \
 public:                                                                       \
  inline auto name(const T &new_##name) -> decltype(*this) { /* NOLINT */      \
    this->name##_ = new_##name;                                                \
    return *this;                                                              \
  }                                                                            \
  inline auto name(T &&new_##name) -> decltype(*this) { /* NOLINT */           \
    this->name##_ = std::move(new_##name);                                     \
    return *this;                                                              \
  }                                                                            \
  DISPATCH_MACRO                                                               \
  inline const T &name() const noexcept { /* NOLINT */ return this->name##_; } \
  DISPATCH_MACRO                                                               \
  inline T &name() noexcept { /* NOLINT */ return this->name##_; }             \
                                                                               \
 private:                                                                      \
  T name##_ /* NOLINT */
