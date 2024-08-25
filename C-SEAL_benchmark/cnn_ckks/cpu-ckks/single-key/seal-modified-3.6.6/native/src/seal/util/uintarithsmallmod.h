// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT license.

#pragma once

#include "seal/modulus.h"
#include "seal/util/defines.h"
#include "seal/util/numth.h"
#include "seal/util/pointer.h"
#include "seal/util/uintarith.h"
#include <cstdint>
#include <type_traits>
using namespace std;

namespace seal
{
    namespace util
    {
        /**
        Returns (operand++) mod modulus.
        Correctness: operand must be at most (2 * modulus -2) for correctness.
        */
        SEAL_NODISCARD inline std::uint64_t increment_uint_mod(std::uint64_t operand, const Modulus &modulus)
        {
#ifdef SEAL_DEBUG
            if (modulus.is_zero())
            {
                throw std::invalid_argument("modulus");
            }
            if (operand > (modulus.value() - 1) << 1)
            {
                throw std::out_of_range("operand");
            }
#endif
            operand++;
            return operand - (modulus.value() &
                              static_cast<std::uint64_t>(-static_cast<std::int64_t>(operand >= modulus.value())));
        }

        /**
        Returns (operand--) mod modulus.
        @param[in] operand Must be at most (modulus - 1).
        */
        SEAL_NODISCARD inline std::uint64_t decrement_uint_mod(std::uint64_t operand, const Modulus &modulus)
        {
#ifdef SEAL_DEBUG
            if (modulus.is_zero())
            {
                throw std::invalid_argument("modulus");
            }
            if (operand >= modulus.value())
            {
                throw std::out_of_range("operand");
            }
#endif
            std::int64_t carry = static_cast<std::int64_t>(operand == 0);
            return operand - 1 + (modulus.value() & static_cast<std::uint64_t>(-carry));
        }

        /**
        Returns (-operand) mod modulus.
        Correctness: operand must be at most modulus for correctness.
        */
        SEAL_NODISCARD inline std::uint64_t negate_uint_mod(std::uint64_t operand, const Modulus &modulus)
        {
#ifdef SEAL_DEBUG
            if (modulus.is_zero())
            {
                throw std::invalid_argument("modulus");
            }
            if (operand >= modulus.value())
            {
                throw std::out_of_range("operand");
            }
#endif
            std::int64_t non_zero = static_cast<std::int64_t>(operand != 0);
            return (modulus.value() - operand) & static_cast<std::uint64_t>(-non_zero);
        }

        /**
        Returns (operand * inv(2)) mod modulus.
        Correctness: operand must be even and at most (2 * modulus - 2) or odd and at most (modulus - 2).
        @param[in] operand Should be at most (modulus - 1).
        */
        SEAL_NODISCARD inline std::uint64_t div2_uint_mod(std::uint64_t operand, const Modulus &modulus)
        {
#ifdef SEAL_DEBUG
            if (modulus.is_zero())
            {
                throw std::invalid_argument("modulus");
            }
            if (operand >= modulus.value())
            {
                throw std::out_of_range("operand");
            }
#endif
            if (operand & 1)
            {
                unsigned long long temp;
                unsigned char carry = add_uint64(operand, modulus.value(), 0, &temp);
                operand = temp >> 1;
                if (carry)
                {
                    return operand | (std::uint64_t(1) << (bits_per_uint64 - 1));
                }
                return operand;
            }
            return operand >> 1;
        }

        /**
        Returns (operand1 + operand2) mod modulus.
        Correctness: (operand1 + operand2) must be at most (2 * modulus - 1).
        */
        SEAL_NODISCARD inline std::uint64_t add_uint_mod(
            std::uint64_t operand1, std::uint64_t operand2, const Modulus &modulus)
        {
#ifdef SEAL_DEBUG
            if (modulus.is_zero())
            {
                throw std::invalid_argument("modulus");
            }
            if (operand1 + operand2 >= modulus.value() << 1)
            {
                throw std::out_of_range("operands");
            }
#endif
            // Sum of operands modulo Modulus can never wrap around 2^64
            operand1 += operand2;
            return SEAL_COND_SELECT(operand1 >= modulus.value(), operand1 - modulus.value(), operand1);
        }

        /**
        Returns (operand1 - operand2) mod modulus.
        Correctness: (operand1 - operand2) must be at most (modulus - 1) and at least (-modulus).
        @param[in] operand1 Should be at most (modulus - 1).
        @param[in] operand2 Should be at most (modulus - 1).
        */
        SEAL_NODISCARD inline std::uint64_t sub_uint_mod(
            std::uint64_t operand1, std::uint64_t operand2, const Modulus &modulus)
        {
#ifdef SEAL_DEBUG
            if (modulus.is_zero())
            {
                throw std::invalid_argument("modulus");
            }

            if (operand1 >= modulus.value())
            {
                throw std::out_of_range("operand1");
            }
            if (operand2 >= modulus.value())
            {
                throw std::out_of_range("operand2");
            }
#endif
            unsigned long long temp;
            std::int64_t borrow = static_cast<std::int64_t>(SEAL_SUB_BORROW_UINT64(operand1, operand2, 0, &temp));
            return static_cast<std::uint64_t>(temp) + (modulus.value() & static_cast<std::uint64_t>(-borrow));
        }

        /**
        Returns input mod modulus. This is not standard Barrett reduction.
        Correctness: modulus must be at most 63-bit.
        @param[in] input Should be at most 128-bit.
        */
        template <typename T, typename = std::enable_if_t<is_uint64_v<T>>>
        SEAL_NODISCARD inline std::uint64_t barrett_reduce_128(const T *input, const Modulus &modulus)
        {
#ifdef SEAL_DEBUG
            if (!input)
            {
                throw std::invalid_argument("input");
            }
            if (modulus.is_zero())
            {
                throw std::invalid_argument("modulus");
            }
#endif
            unsigned long long tmp1, tmp2[2], tmp3, carry;
            const std::uint64_t *const_ratio_ori = modulus.const_ratio().data();

            // Multiply input and const_ratio
            // Round 1
            multiply_uint64_hw64(input[0], const_ratio_ori[0], &carry);

            multiply_uint64(input[0], const_ratio_ori[1], tmp2);
            tmp3 = tmp2[1] + add_uint64(tmp2[0], carry, &tmp1);

            // Round 2
            multiply_uint64(input[1], const_ratio_ori[0], tmp2);
            carry = tmp2[1] + add_uint64(tmp1, tmp2[0], &tmp1);

            // This is all we care about
            tmp1 = input[1] * const_ratio_ori[1] + tmp3 + carry;

            // Barrett subtraction
            tmp3 = input[0] - tmp1 * modulus.value();

            // One more subtraction is enough
            uint64_t result_origin = SEAL_COND_SELECT(tmp3 >= modulus.value(), tmp3 - modulus.value(), tmp3);

            return result_origin;
        }

        /**
        Returns input mod modulus. This is not standard Barrett reduction.
        Correctness: modulus must be at most 63-bit.
        @param[in] input Should be at most 128-bit.
        */
        template <typename T, typename = std::enable_if_t<is_uint64_v<T>>>
        SEAL_NODISCARD inline std::uint64_t barrett_reduce_128_nonlazy(const T *input, const Modulus &modulus)
        {
#ifdef SEAL_DEBUG
            if (!input)
            {
                throw std::invalid_argument("input");
            }
            if (modulus.is_zero())
            {
                throw std::invalid_argument("modulus");
            }
#endif
            // unsigned long long tmp1, tmp2[2], tmp3, carry;
            // const std::uint64_t *const_ratio_ori = modulus.const_ratio().data();

            // // Multiply input and const_ratio
            // // Round 1
            // multiply_uint64_hw64(input[0], const_ratio_ori[0], &carry);

            // multiply_uint64(input[0], const_ratio_ori[1], tmp2);
            // tmp3 = tmp2[1] + add_uint64(tmp2[0], carry, &tmp1);

            // // Round 2
            // multiply_uint64(input[1], const_ratio_ori[0], tmp2);
            // carry = tmp2[1] + add_uint64(tmp1, tmp2[0], &tmp1);

            // // This is all we care about
            // tmp1 = input[1] * const_ratio_ori[1] + tmp3 + carry;

            // // Barrett subtraction
            // tmp3 = input[0] - tmp1 * modulus.value();

            // // One more subtraction is enough
            // uint64_t result_origin = SEAL_COND_SELECT(tmp3 >= modulus.value(), tmp3 - modulus.value(), tmp3);

            // Reduces input using base 2^64 Barrett reduction  Chao Added
            // input allocation size must be 128 bits
            unsigned long long quo;
            const std::uint64_t *const_ratio = modulus.const_ratio128().data();

            std::uint64_t TH[2]{ 0, 0 };
            // TH[0] = input[0];
            // TH[1] = input[1];
            const uint64_t *input_as_uint64 = reinterpret_cast<const uint64_t*>(input);
            // Right shift 'input' by (modulus bit count - 2) bits and store result in 'TH'
            right_shift_uint128(input_as_uint64, modulus.bit_count() - 2, TH);
            // TH * M >> N+3
            multiply_uint64_hw64(TH[0], const_ratio[0], &quo);
            uint64_t tmp = quo * modulus.value();
            // Calculate the candidate for modulus
            uint64_t Z = input[0] - tmp;
            // Modular correction step to ensure the result is in range [0, modulus)
            uint64_t result = (Z >= modulus.value()) ? Z - modulus.value() : Z;
            return result;
            // bool erro_occur = false;
            // if (result == result_origin){
            //     return result;
            // }
            // else {
            //     erro_occur = true;

            //     int bit_size = get_significant_bit_count(result_origin);
            //     std::cout << "origin result: "<<std::hex<<result_origin<<std::dec<<", bitsize: "<<bit_size<<std::endl;
            //     bit_size = get_significant_bit_count(result);
            //     std::cout << "new result: "<<std::hex<<result<<std::dec<<", bitsize: "<<bit_size<<std::endl;
            //     cout << "input[0] is lower part, input[0]: " <<hex<<input[0]<<", input[1]: "<<input[1]<<endl;
            //     cout << "After right shift TH[0]: "<<TH[0]<<", TH[1]: "<<TH[1]<<endl;
            //     cout << "After multiply_uint64_hw64, TH[0]: "<<TH[0]<<", const_ratio[0]: "<<const_ratio[0]<<", quo: "<<quo<<endl;
            //     cout << "quo: "<<quo<<"* modulus: "<<modulus.value()<<" = tmp: "<<tmp<<endl;
            //     cout << "Z :"<<Z<<endl;

            //     if (erro_occur){ throw std::invalid_argument("invalid modulu reduction");  }
            //     return result_origin;
            // }
        }



        /**
        Returns input mod modulus. This is not standard Barrett reduction.
        Correctness: modulus must be at most 63-bit.
        */
        template <typename T, typename = std::enable_if_t<is_uint64_v<T>>>
        SEAL_NODISCARD inline std::uint64_t barrett_reduce_64(T input, const Modulus &modulus)
        {
#ifdef SEAL_DEBUG
            if (modulus.is_zero())
            {
                throw std::invalid_argument("modulus");
            }
#endif
            // Reduces input using base 2^64 Barrett reduction
            // floor(2^64 / mod) == floor( floor(2^128 / mod) )
            unsigned long long tmp[2];
            const std::uint64_t *const_ratio = modulus.const_ratio().data();
            multiply_uint64_hw64(input, const_ratio[1], tmp + 1);
            // uint64_t lo64 = tmp[1] * modulus.value();
            // Barrett subtraction
            tmp[0] = input - tmp[1] * modulus.value();
            // tmp[0] = input - lo64;

            // One more subtraction is enough
            return SEAL_COND_SELECT(tmp[0] >= modulus.value(), tmp[0] - modulus.value(), tmp[0]);
        }

        /**
        Returns (operand1 * operand2) mod modulus.
        Correctness: Follows the condition of barrett_reduce_128.
        */
       //Chao added, input is not limited to [0,M), may be a lazy reduction
        SEAL_NODISCARD inline std::uint64_t multiply_uint_mod(
            std::uint64_t operand1, std::uint64_t operand2, const Modulus &modulus)
        {
#ifdef SEAL_DEBUG
            if (modulus.is_zero())
            {
                throw std::invalid_argument("modulus");
            }
#endif
            unsigned long long z[2];
            multiply_uint64(operand1, operand2, z);
            return barrett_reduce_128(z, modulus); // right
        }

        SEAL_NODISCARD inline std::uint64_t multiply_uint_mod_nonlazy(
            std::uint64_t operand1, std::uint64_t operand2, const Modulus &modulus)
        {
#ifdef SEAL_DEBUG
            if (modulus.is_zero())
            {
                throw std::invalid_argument("modulus");
            }
#endif
            unsigned long long z[2];
            multiply_uint64(operand1, operand2, z);
            return barrett_reduce_128_nonlazy(z, modulus); // right
        }

        /**
        This struct contains a operand and a precomputed quotient: (operand << 64) / modulus, for a specific modulus.
        When passed to multiply_uint_mod, a faster variant of Barrett reduction will be performed.
        Operand must be less than modulus.
        */
        struct MultiplyUIntModOperand
        {
            std::uint64_t operand;
            std::uint64_t quotient;

            void set_quotient(const Modulus &modulus)
            {
#ifdef SEAL_DEBUG
                if (operand >= modulus.value())
                {
                    throw std::invalid_argument("input must be less than modulus");
                }
#endif
                std::uint64_t wide_quotient[2]{ 0, 0 };
                std::uint64_t wide_coeff[2]{ 0, operand };
                divide_uint128_inplace(wide_coeff, modulus.value(), wide_quotient);
                quotient = wide_quotient[0];
            }

            void set(std::uint64_t new_operand, const Modulus &modulus)
            {
#ifdef SEAL_DEBUG
                if (new_operand >= modulus.value())
                {
                    throw std::invalid_argument("input must be less than modulus");
                }
#endif
                operand = new_operand;
                set_quotient(modulus);
            }
        };

        /**
        Returns x * y mod modulus.
        This is a highly-optimized variant of Barrett reduction.
        Correctness: modulus should be at most 63-bit, and y must be less than modulus.
        */
        SEAL_NODISCARD inline std::uint64_t multiply_uint_mod(
            std::uint64_t x, MultiplyUIntModOperand y, const Modulus &modulus)
        {
#ifdef SEAL_DEBUG
            if (y.operand >= modulus.value())
            {
                throw std::invalid_argument("operand y must be less than modulus");
            }
#endif
            unsigned long long tmp1, tmp2;
            const std::uint64_t p = modulus.value();
            multiply_uint64_hw64(x, y.quotient, &tmp1);
            tmp2 = y.operand * x - tmp1 * p;
            return SEAL_COND_SELECT(tmp2 >= p, tmp2 - p, tmp2);
        }

        /**
        Returns x * y mod modulus or x * y mod modulus + modulus.
        This is a highly-optimized variant of Barrett reduction and reduce to [0, 2 * modulus - 1].
        Correctness: modulus should be at most 63-bit, and y must be less than modulus.
        */
        SEAL_NODISCARD inline std::uint64_t multiply_uint_mod_lazy(
            std::uint64_t x, MultiplyUIntModOperand y, const Modulus &modulus)
        {
#ifdef SEAL_DEBUG
            if (y.operand >= modulus.value())
            {
                throw std::invalid_argument("operand y must be less than modulus");
            }
#endif
            unsigned long long tmp1;
            const std::uint64_t p = modulus.value();
            multiply_uint64_hw64(x, y.quotient, &tmp1);
            return y.operand * x - tmp1 * p;
        }

        /**
        Returns value[0] = value mod modulus.
        Correctness: Follows the condition of barrett_reduce_128.
        */
        inline void modulo_uint_inplace(std::uint64_t *value, std::size_t value_uint64_count, const Modulus &modulus)
        {
#ifdef SEAL_DEBUG
            if (!value)
            {
                throw std::invalid_argument("value");
            }
            if (!value_uint64_count)
            {
                throw std::invalid_argument("value_uint64_count");
            }
#endif

            if (value_uint64_count == 1)
            {
                if (*value < modulus.value())
                {
                    return;
                }
                else
                {
                    *value = barrett_reduce_64(*value, modulus);
                }
            }

            // Starting from the top, reduce always 128-bit blocks
            for (std::size_t i = value_uint64_count - 1; i--;)
            {
                value[i] = barrett_reduce_128_nonlazy(value + i, modulus); // added
                value[i + 1] = 0;
            }
        }

        /**
        Returns value mod modulus.
        Correctness: Follows the condition of barrett_reduce_128.
        */
        SEAL_NODISCARD inline std::uint64_t modulo_uint(
            const std::uint64_t *value, std::size_t value_uint64_count, const Modulus &modulus)
        {
#ifdef SEAL_DEBUG
            if (!value && value_uint64_count)
            {
                throw std::invalid_argument("value");
            }
            if (!value_uint64_count)
            {
                throw std::invalid_argument("value_uint64_count");
            }
#endif
            if (value_uint64_count == 1)
            {
                // If value < modulus no operation is needed
                if (*value < modulus.value())
                    return *value;
                else
                    return barrett_reduce_64(*value, modulus);
            }

            // Temporary space for 128-bit reductions
            uint64_t temp[2]{ 0, value[value_uint64_count - 1] };
            for (size_t k = value_uint64_count - 1; k--;)
            {
                temp[0] = value[k];
                temp[1] = barrett_reduce_128(temp, modulus); // 128bit input of barret reduction
            }

            // Save the result modulo i-th prime
            return temp[1];
        }

        /**
        Returns (operand1 * operand2) + operand3 mod modulus.
        Correctness: Follows the condition of barrett_reduce_128.
        */
        inline std::uint64_t multiply_add_uint_mod(
            std::uint64_t operand1, std::uint64_t operand2, std::uint64_t operand3, const Modulus &modulus)
        {
            // Lazy reduction
            unsigned long long temp[2];
            multiply_uint64(operand1, operand2, temp);
            temp[1] += add_uint64(temp[0], operand3, temp);
            return barrett_reduce_128_nonlazy(temp, modulus);//third change
        }

        /**
        Returns (operand1 * operand2) + operand3 mod modulus.
        Correctness: Follows the condition of multiply_uint_mod.
        */
        inline std::uint64_t multiply_add_uint_mod(
            std::uint64_t operand1, MultiplyUIntModOperand operand2, std::uint64_t operand3, const Modulus &modulus)
        {
            return add_uint_mod(
                multiply_uint_mod(operand1, operand2, modulus), barrett_reduce_64(operand3, modulus), modulus);
        }

        inline bool try_invert_uint_mod(std::uint64_t operand, const Modulus &modulus, std::uint64_t &result)
        {
            return try_invert_uint_mod(operand, modulus.value(), result);
        }

        /**
        Returns operand^exponent mod modulus.
        Correctness: Follows the condition of barrett_reduce_128.
        */
        SEAL_NODISCARD std::uint64_t exponentiate_uint_mod(
            std::uint64_t operand, std::uint64_t exponent, const Modulus &modulus);

        /**
        Computes numerator = numerator mod modulus, quotient = numerator / modulus.
        Correctness: Follows the condition of barrett_reduce_128.
        */
        void divide_uint_mod_inplace(
            std::uint64_t *numerator, const Modulus &modulus, std::size_t uint64_count, std::uint64_t *quotient,
            MemoryPool &pool);

        /**
        Computes <operand1, operand2> mod modulus.
        Correctness: Follows the condition of barrett_reduce_128.
        */
        SEAL_NODISCARD std::uint64_t dot_product_mod(
            const std::uint64_t *operand1, const std::uint64_t *operand2, std::size_t count, const Modulus &modulus);
    } // namespace util
} // namespace seal
