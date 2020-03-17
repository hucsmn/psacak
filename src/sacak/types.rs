use std::mem::size_of;

use std::fmt::{Debug, Display};
use std::num::{NonZeroU128, NonZeroU16, NonZeroU32, NonZeroU64, NonZeroU8, NonZeroUsize};
use std::ops::{Add, BitAnd, BitOr, BitXor, Not, Shl, Shr, Sub};
use std::ops::{AddAssign, BitAndAssign, BitOrAssign, BitXorAssign, ShlAssign, ShrAssign, SubAssign};

// Cheap unsigned integer cast.
pub trait As<T: Copy>: Copy {
    fn r#as(self) -> T;
}

macro_rules! impl_as {
    ($($t1:ty => $( $t2:ty ),* ; )*) => {
        $( $(
            impl As<$t2> for $t1 {
                #[inline(always)]
                fn r#as(self) -> $t2 {
                    self as $t2
                }
            }
        )* )*
    };
}

impl_as! {
    u8 => u8, u16, u32, u64, u128, usize;
    u16 => u8, u16, u32, u64, u128, usize;
    u32 => u8, u16, u32, u64, u128, usize;
    u64 => u8, u16, u32, u64, u128, usize;
    u128 => u8, u16, u32, u64, u128, usize;
    usize => u8, u16, u32, u64, u128, usize;
}

/// Unsigned integers with basic arithmetic operations.
pub trait Uint:
    Copy
    + Eq
    + Ord
    + Debug
    + Display
    + ToString
    + As<u8>
    + As<u16>
    + As<u32>
    + As<u64>
    + As<u128>
    + As<usize>
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

    type NonZero: NonZeroUint;

    fn into_nonzero(self) -> Option<Self::NonZero>;

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
    ($($uint:ident($nonzero:ident);)*) => {
        $(
            impl Uint for $uint {
                const ZERO: Self = 0;
                const ONE: Self = 1;
                const MAX: Self = std::$uint::MAX;
                const SIZE: usize = size_of::<$uint>();
                const BIT_WIDTH: u8 = 8 * (size_of::<$uint>() as u8);

                type NonZero = $nonzero;

                #[inline(always)]
                fn into_nonzero(self) -> Option<Self::NonZero> {
                    Self::NonZero::new(self)
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

impl_uint! {
    u8(NonZeroU8);
    u16(NonZeroU16);
    u32(NonZeroU32);
    u64(NonZeroU64);
    u128(NonZeroU128);
    usize(NonZeroUsize);
}

/// Non-zero unsigned integer type.
pub trait NonZeroUint: Copy + Eq + Ord {
    type Zeroable: Uint;

    fn from_zeroable(x: Self::Zeroable) -> Option<Self>;
    fn into_zeroable(self) -> Self::Zeroable;
}

macro_rules! impl_nonzero {
    ($($nonzero:ident($uint:ident);)*) => {
        $(
            impl NonZeroUint for $nonzero {
                type Zeroable = $uint;

                #[inline(always)]
                fn from_zeroable(x: $uint) -> Option<Self> {
                    Self::new(x)
                }

                #[inline(always)]
                fn into_zeroable(self) -> $uint {
                    self.get()
                }
            }
        )*
    };
}

impl_nonzero! {
    NonZeroU8(u8);
    NonZeroU16(u16);
    NonZeroU32(u32);
    NonZeroU64(u64);
    NonZeroU128(u128);
    NonZeroUsize(usize);
}

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
