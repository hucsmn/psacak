use std::mem::size_of;

use std::fmt::{Debug, Display};
use std::ops::{Add, BitAnd, BitOr, BitXor, Not, Shl, Shr, Sub};
use std::ops::{
    AddAssign, BitAndAssign, BitOrAssign, BitXorAssign, ShlAssign, ShrAssign, SubAssign,
};

/// Unsigned integers with basic arithmetic operations.
pub trait Uint:
    Copy
    + Eq
    + Ord
    + Debug
    + Display
    + ToString
    + Add<Self, Output = Self>
    + Sub<Self, Output = Self>
    + Not<Output = Self>
    + BitAnd<Self, Output = Self>
    + BitOr<Self, Output = Self>
    + BitXor<Self, Output = Self>
    + Shl<u8, Output = Self>
    + Shr<u8, Output = Self>
    + AddAssign<Self>
    + SubAssign<Self>
    + BitAndAssign<Self>
    + BitOrAssign<Self>
    + BitXorAssign<Self>
    + ShlAssign<u8>
    + ShrAssign<u8>
{
    const ZERO: Self;
    const ONE: Self;
    const MAX: Self;
    const SIZE: usize;
    const BIT_WIDTH: u8;

    fn from_u64(x: u64) -> Self;
    fn as_u64(self) -> u64;

    fn wrapping_add(self, other: Self) -> Self;
    fn saturating_add(self, other: Self) -> Self;
    fn checked_add(self, other: Self) -> Option<Self>;

    fn wrapping_sub(self, other: Self) -> Self;
    fn saturating_sub(self, other: Self) -> Self;
    fn checked_sub(self, other: Self) -> Option<Self>;
}

macro_rules! forward_binops_for_impl_uint {
    ($uint:ident => $($method:ident -> $ret:ty),*) => {
        $(
            #[inline(always)]
            fn $method(self, other: Self) -> $ret {
                <$uint>::$method(self, other)
            }
        )*
    };
}

macro_rules! impl_uint {
    ($($uint:ident),*) => {
        $(
            impl Uint for $uint {
                const ZERO: Self = 0;
                const ONE: Self = 1;
                const MAX: Self = std::$uint::MAX;
                const SIZE: usize = size_of::<$uint>();
                const BIT_WIDTH: u8 = 8 * (size_of::<$uint>() as u8);

                #[inline(always)]
                fn from_u64(x: u64) -> Self {
                    x as $uint
                }

                #[inline(always)]
                fn as_u64(self) -> u64 {
                    self as u64
                }

                forward_binops_for_impl_uint! { $uint =>
                    wrapping_add -> Self,
                    saturating_add -> Self,
                    checked_add -> Option<Self>,
                    wrapping_sub -> Self,
                    saturating_sub -> Self,
                    checked_sub -> Option<Self>
                }
            }
        )*
    };
}

impl_uint!(u8, u32, usize);

/// Types that could be casted into usize.
pub trait AsIndex: Copy {
    fn as_index(self) -> usize;
}

macro_rules! impl_as_index {
    ($($uint:ty),*) => {
        $(
            impl AsIndex for $uint {
                #[inline(always)]
                fn as_index(self) -> usize {
                    self as usize
                }
            }
        )*
    };
}

impl_as_index!(u8, u32, usize);

/// Text character type.
pub trait SacaChar: Uint + AsIndex {}

impl SacaChar for u8 {}

impl SacaChar for u32 {}

/// Suffix array index type.
pub trait SacaIndex: Uint + AsIndex {
    fn from_index(idx: usize) -> Self;
}

macro_rules! impl_saca_index {
    ($($uint:ty),*) => {
        $(
            impl SacaIndex for $uint {
                #[inline(always)]
                fn from_index(idx: usize) -> Self {
                    idx as $uint
                }
            }
        )*
    };
}

impl_saca_index!(u32);
