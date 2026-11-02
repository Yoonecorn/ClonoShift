// SPDX-License-Identifier: GPL-3.0-or-later
//! Opt-in stage wall durations on stderr, outside all scientific artifacts.
use std::time::Instant;
pub(crate) struct Profile {
    previous: Option<Instant>,
    workflow: &'static str,
}
impl Profile {
    pub(crate) fn new(workflow: &'static str) -> Self {
        Self {
            previous: (std::env::var_os("CLONOSHIFT_PROFILE").as_deref()
                == Some(std::ffi::OsStr::new("1")))
            .then(Instant::now),
            workflow,
        }
    }
    pub(crate) fn mark(&mut self, stage: &'static str) {
        if let Some(previous) = self.previous.as_mut() {
            let now = Instant::now();
            eprintln!(
                "clonoshift-profile {} {} {:.9}",
                self.workflow,
                stage,
                now.duration_since(*previous).as_secs_f64()
            );
            *previous = now;
        }
    }
}
