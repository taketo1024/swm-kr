// swift-tools-version:5.3
// The swift-tools-version declares the minimum version of Swift required to build this package.

import PackageDescription

let package = Package(
    name: "swm-kr",
    products: [
        .library(
            name: "SwmKR",
            targets: ["SwmKR"]
        ),
    ],
    dependencies: [
        .package(
			url: "https://github.com/taketo1024/swm-core.git",
			from:"1.3.1"
//            path: "../swm-core/"
		),
        .package(
            url: "https://github.com/taketo1024/swm-khovanov.git",
            from: "1.2.0"
//            path: "../swm-khovanov/"
        ),
    ],
    targets: [
        .target(
            name: "SwmKR",
            dependencies: [
                .product(name: "SwmCore", package: "swm-core"),
                .product(name: "SwmKhovanov", package: "swm-khovanov"),
            ]
        ),
        .testTarget(
            name: "SwmKRTests",
            dependencies: ["SwmKR"]
		),
    ]
)
