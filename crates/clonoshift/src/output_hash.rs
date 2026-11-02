// SPDX-License-Identifier: GPL-3.0-or-later
//! Bounded streaming SHA-256. Platform failures refuse output publication.
use sha2::{Digest, Sha256};
use std::io;

pub(crate) struct OutputHasher {
    portable: Sha256,
    #[cfg(all(target_os = "macos", target_arch = "aarch64"))]
    native: Option<apple::Context>,
}
impl OutputHasher {
    pub(crate) fn new() -> Self {
        Self {
            portable: Sha256::new(),
            #[cfg(all(target_os = "macos", target_arch = "aarch64"))]
            native: apple::Context::new(),
        }
    }
    pub(crate) fn update(&mut self, bytes: &[u8]) -> io::Result<()> {
        #[cfg(all(target_os = "macos", target_arch = "aarch64"))]
        if let Some(native) = self.native.as_mut() {
            return native.update(bytes);
        }
        self.portable.update(bytes);
        Ok(())
    }
    pub(crate) fn finish(self) -> io::Result<String> {
        #[cfg(all(target_os = "macos", target_arch = "aarch64"))]
        if let Some(native) = self.native {
            let bytes = native.finish()?;
            return Ok(bytes.iter().map(|b| format!("{b:02x}")).collect());
        }
        Ok(format!("{:x}", self.portable.finalize()))
    }
}

#[cfg(all(target_os = "macos", target_arch = "aarch64"))]
mod apple {
    use std::io;
    // Public CommonDigest.h ABI: CC_LONG=uint32_t; CC_SHA256_CTX contains
    // count[2], hash[8], wbuf[16]. No Apple implementation is copied.
    #[repr(C)]
    pub(super) struct Context {
        count: [u32; 2],
        hash: [u32; 8],
        buffer: [u32; 16],
    }
    #[link(name = "System")]
    unsafe extern "C" {
        fn CC_SHA256_Init(context: *mut Context) -> libc::c_int;
        fn CC_SHA256_Update(
            context: *mut Context,
            bytes: *const libc::c_void,
            length: u32,
        ) -> libc::c_int;
        fn CC_SHA256_Final(digest: *mut u8, context: *mut Context) -> libc::c_int;
    }
    impl Context {
        pub(super) fn new() -> Option<Self> {
            let mut context = Self {
                count: [0; 2],
                hash: [0; 8],
                buffer: [0; 16],
            };
            // SAFETY: The correctly laid out owned context is fully writable.
            (unsafe { CC_SHA256_Init(&mut context) } == 1).then_some(context)
        }
        pub(super) fn update(&mut self, bytes: &[u8]) -> io::Result<()> {
            for chunk in bytes.chunks(u32::MAX as usize) {
                // SAFETY: Context was initialized, chunk remains live, length
                // is checked by chunking, and the API retains no input pointer.
                if unsafe { CC_SHA256_Update(self, chunk.as_ptr().cast(), chunk.len() as u32) } != 1
                {
                    return Err(io::Error::other("Apple SHA-256 update failed"));
                }
            }
            Ok(())
        }
        pub(super) fn finish(mut self) -> io::Result<[u8; 32]> {
            let mut digest = [0; 32];
            // SAFETY: Initialized owned context and writable32byte output;
            // neither pointer is retained or shared with another hasher.
            if unsafe { CC_SHA256_Final(digest.as_mut_ptr(), &mut self) } != 1 {
                return Err(io::Error::other("Apple SHA-256 finalization failed"));
            }
            Ok(digest)
        }
    }
    #[cfg(test)]
    #[test]
    fn sdk_context_layout_matches_c_abi() {
        assert_eq!(std::mem::size_of::<Context>(), 104);
        assert_eq!(std::mem::align_of::<Context>(), 4);
        assert_eq!(std::mem::offset_of!(Context, hash), 8);
        assert_eq!(std::mem::offset_of!(Context, buffer), 40);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn streaming_matches_independent_portable_hash_across_chunk_boundaries() {
        #[cfg(all(target_os = "macos", target_arch = "aarch64"))]
        assert!(
            OutputHasher::new().native.is_some(),
            "native API must be exercised on this host"
        );
        for length in [0, 1, 55, 56, 63, 64, 65, 129, 1_048_577] {
            let input: Vec<u8> = (0..length).map(|i| (i * 73 + 17) as u8).collect();
            let expected = format!("{:x}", Sha256::digest(&input));
            for chunk in [1, 7, 55, 64, 127, 131_072] {
                let mut hasher = OutputHasher::new();
                for part in input.chunks(chunk) {
                    hasher.update(part).unwrap();
                }
                assert_eq!(hasher.finish().unwrap(), expected);
            }
        }
        #[cfg(all(target_os = "macos", target_arch = "aarch64"))]
        assert_eq!(std::mem::size_of::<apple::Context>(), 104);
    }
}
