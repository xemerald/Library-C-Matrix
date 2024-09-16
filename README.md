# C Library for Matrix calculation

## Dependencies
Nothing special, it can be compiled under standard C library.

## Supported Platforms
- Linux
- MacOS (**Not pass testing yet**)

## Building & Installation
- Linux/MacOS
	- Simply run `make` or `make default` to get the static library `libmatrix.a` with usage of CPU intrinsics
	- Or run `make static` to get the **naive** library `libmatrix.a` without usage of CPU intrinsics
	- Or run `make share` or `make dynamic` to get the **naive** shared library `libmatrix.so` without usage of CPU intrinsics
	- Then you can run `sudo make install` to install this library to `/usr/local`

## Unit Testing
- Linux/MacOS
	- Simply run `make test` to run the unit test

## Configuration

## License
Copyright&copy; 2023-2024 Benjamin Ming Yang

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at

[Apache License 2.0](http://www.apache.org/licenses/LICENSE-2.0)

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.
