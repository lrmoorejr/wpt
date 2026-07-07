# Copyright 2026 L. Richard Moore Jr.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

cmake_minimum_required(VERSION 3.18)
project(wpt VERSION 1.0.0 LANGUAGES CXX)

if(NOT CMAKE_CXX_STANDARD)
	# WPT.hpp's actual floor is std::make_unique (C++14) -- only default to it when
	# nothing higher is already set. Built standalone this verifies that floor via
	# the assert()-fallback path (no Ensure.hpp on the include path); embedded in a
	# C++20-or-newer parent build, it inherits that higher standard instead, so the
	# real Ensure.hpp (which itself requires C++20) can be used when found.
	set(CMAKE_CXX_STANDARD 14)
endif()
# Catches member-initializer-list order mismatching declaration order (C++
# initializes in declaration order regardless of initializer-list order) --
# the root cause of a past intermittent segfault in this header.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wreorder -Wall -Wextra -Wpedantic -Werror")

find_package(Catch2 3 REQUIRED)

add_executable(WPT-test WPT-test.cpp WPT.hpp)
target_link_libraries(WPT-test PRIVATE Catch2::Catch2WithMain)

enable_testing()
include(Catch)
catch_discover_tests(WPT-test)
