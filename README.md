# Purpose
We have publicly forked the Zcash library beacuase we are refactoring bellman. In its current state, the library is difficult to understand if you are not familiar with both the original Groth16 paper, as well as the work by the Zcash and Sonic team on extending the capabilities of the implementation. One consequence of this is that vetting the library for bugs and errors is considerably harder, because we are attempting to understand 200+ line functions with ambiguous names deeply connected to the academic works. Maintaining the work is also somewhat harder, because the original team is required to understand particular design decisions.

The most up-to-date work is in the 'refactor' branch, and we also have a 'better test' fork which breaks up the original tests for better readability. 

Feel free to raise issues if you have an improvement (however small) or notice that we have made an error, especially because there are parts of the library where we have had to make assumptions because the original work has elided or not documented components of the design.

## Zcash Rust crates

This repository contains a fork of a (work-in-progress) set of Rust crates for
working with Zcash.

### Security Warnings

These libraries are currently under development and have not been fully-reviewed. 

### License

All code in this workspace is licensed under either of

 * Apache License, Version 2.0, ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
 * MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.

#### Contribution

Unless you explicitly state otherwise, any contribution intentionally
submitted for inclusion in the work by you, as defined in the Apache-2.0
license, shall be dual licensed as above, without any additional terms or
conditions.
