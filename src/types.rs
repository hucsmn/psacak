use std::mem::size_of;

use std::fmt::{Debug, Display};
use std::ops::{Add, BitAnd, BitOr, BitXor, Not, Shl, Shr, Sub};
use std::ops::{AddAssign, BitAndAssign, BitOrAssign, BitXorAssign, ShlAssign, ShrAssign, SubAssign};
use std::sync::atomic::Ordering;

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
    + Sync
    + Send
    + Eq
    + Ord
    + Default
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
    const HIGHEST_BIT: Self;
    const LOWER_BITS: Self;

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
                const HIGHEST_BIT: Self = 1 << (Self::BIT_WIDTH - 1);
                const LOWER_BITS: Self = !Self::HIGHEST_BIT;

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

impl_uint!(u8, u16, u32, u64, u128, usize);

/// Mark unsigned integer type has atomic version.
pub trait HasAtomic: Uint {
    type Atomic: Atomic<Self> + From<Self>;
}

/// Atomic unsigned integers.
pub trait Atomic<T: Uint>: Sync + Send + Debug {
    fn new(val: T) -> Self;
    fn get_mut(&mut self) -> &mut T;
    fn into_inner(self) -> T;

    fn load(&self, order: Ordering) -> T;
    fn store(&self, val: T, order: Ordering);
    fn swap(&self, val: T, order: Ordering) -> T;
    fn compare_and_swap(&self, current: T, new: T, order: Ordering) -> T;
    fn compare_exchange(&self, current: T, new: T, success: Ordering, failure: Ordering) -> Result<T, T>;
    fn compare_exchange_weak(&self, current: T, new: T, success: Ordering, failure: Ordering) -> Result<T, T>;

    fn fetch_add(&self, val: T, order: Ordering) -> T;
    fn fetch_sub(&self, val: T, order: Ordering) -> T;
    fn fetch_and(&self, val: T, order: Ordering) -> T;
    fn fetch_nand(&self, val: T, order: Ordering) -> T;
    fn fetch_or(&self, val: T, order: Ordering) -> T;
    fn fetch_xor(&self, val: T, order: Ordering) -> T;
}

macro_rules! forward_binops_for_impl_atomic {
    ($uint:ident $atomic:ident => $($method:ident),*) => {
        $(
            #[inline(always)]
            fn $method(&self, val: $uint, order: Ordering) -> $uint {
                $atomic::$method(self, val, order)
            }
        )*
    };
}

macro_rules! impl_atomic {
    ($($uint:ident $atomic:ident),*) => {
        $(
            impl HasAtomic for $uint {
                type Atomic = $atomic;
            }

            impl Atomic<$uint> for $atomic {
                #[inline(always)]
                fn new(val: $uint) -> Self {
                    $atomic::new(val)
                }
                #[inline(always)]
                fn get_mut(&mut self) -> &mut $uint {
                    $atomic::get_mut(self)
                }
                #[inline(always)]
                fn into_inner(self) -> $uint {
                    $atomic::into_inner(self)
                }

                #[inline(always)]
                fn load(&self, order: Ordering) -> $uint {
                    $atomic::load(self, order)
                }
                #[inline(always)]
                fn store(&self, val: $uint, order: Ordering) {
                    $atomic::store(self, val, order);
                }
                #[inline(always)]
                fn swap(&self, val: $uint, order: Ordering) -> $uint {
                    $atomic::swap(self, val, order)
                }
                #[inline(always)]
                fn compare_and_swap(&self, current: $uint, new: $uint, order: Ordering) -> $uint {
                    $atomic::compare_and_swap(self, current, new, order)
                }
                #[inline(always)]
                fn compare_exchange(&self, current: $uint, new: $uint, success: Ordering, failure: Ordering) -> Result<$uint, $uint> {
                    $atomic::compare_exchange(self, current, new, success, failure)
                }
                #[inline(always)]
                fn compare_exchange_weak(&self, current: $uint, new: $uint, success: Ordering, failure: Ordering) -> Result<$uint, $uint> {
                    $atomic::compare_exchange_weak(self, current, new, success, failure)
                }

                forward_binops_for_impl_atomic!($uint $atomic => fetch_add, fetch_sub, fetch_and, fetch_nand, fetch_or, fetch_xor);
            }
        )*
    };
}

// TODO: rewrite when `#[cfg(target_has_atomic = "<width>")]` comes to stable.
cfg_if! {
    if #[cfg(any(
        target_arch = "x86_64",
        target_arch = "x86",
        target_arch = "aarch64",
        target_arch = "mips64",
        target_arch = "powerpc64",
        target_arch = "riscv64",
        target_arch = "sparc",
        target_arch = "sparc64",
    ))]
    {
        // platforms support atomic u8..u64.
        use std::sync::atomic::{AtomicU8, AtomicU16, AtomicU32, AtomicU64, AtomicUsize};
        impl_atomic!(u8 AtomicU8, u16 AtomicU16, u32 AtomicU32, u64 AtomicU64, usize AtomicUsize);
    } else if #[cfg(any(
        target_arch = "arm",
        target_arch = "mips",
        target_arch = "powerpc",
        target_arch = "wasm32",
    ))]
    {
        // platforms support atomic u8..u32.
        use std::sync::atomic::{AtomicU8, AtomicU16, AtomicU32, AtomicUsize};
        impl_atomic!(u8 AtomicU8, u16 AtomicU16, u32 AtomicU32, usize AtomicUsize);
    }
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

/// Text character type.
pub trait SacaChar: Uint + HasAtomic + AsIndex {}

macro_rules! impl_saca_char {
    ($($uint:ty),*) => {
        $(
            impl SacaChar for $uint {}
        )*
    };
}

/// Suffix array index type.
pub trait SacaIndex: Uint + HasAtomic + AsIndex {
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

cfg_if! {
    if #[cfg(target_pointer_width="64")] {
        impl_as_index!(u8, u16, u32, u64, usize);
        impl_saca_char!(u8, u16, u32, u64);
        impl_saca_index!(u32, u64);
    } else if #[cfg(target_pointer_width="32")] {
        impl_as_index!(u8, u16, u32, usize);
        impl_saca_char!(u8, u16, u32);
        impl_saca_index!(u32);
    }
}
