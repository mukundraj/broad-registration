
export function draw_bounding_box(origin, corner, scene){

    const myLines = [
        [new BABYLON.Vector3(origin[0], origin[1], origin[2]),
            new BABYLON.Vector3(corner[0], origin[1], origin[2])
        ],
        [new BABYLON.Vector3(origin[0], origin[1], origin[2]),
            new BABYLON.Vector3(origin[0], corner[1], origin[2])
        ],
        [new BABYLON.Vector3(origin[0], origin[1], origin[2]),
            new BABYLON.Vector3(origin[0], origin[1], corner[2])
        ],
        [new BABYLON.Vector3(origin[0], corner[1], origin[2]),
            new BABYLON.Vector3(corner[0], corner[1], origin[2])
        ],
        [new BABYLON.Vector3(corner[0], corner[1], origin[2]),
            new BABYLON.Vector3(corner[0], origin[1], origin[2])
        ],
        [new BABYLON.Vector3(origin[0], origin[1], corner[2]),
            new BABYLON.Vector3(origin[0], corner[1], corner[2])
        ],
        [new BABYLON.Vector3(origin[0], corner[1], origin[2]),
            new BABYLON.Vector3(origin[0], corner[1], corner[2])
        ],
        [new BABYLON.Vector3(origin[0], origin[1], corner[2]),
            new BABYLON.Vector3(corner[0], origin[1], corner[2])
        ],
        [new BABYLON.Vector3(corner[0], corner[1], origin[2]),
            new BABYLON.Vector3(corner[0], corner[1], corner[2])
        ],
        [new BABYLON.Vector3(origin[0], corner[1], corner[2]),
            new BABYLON.Vector3(corner[0], corner[1], corner[2])
        ],
        [new BABYLON.Vector3(corner[0], origin[1], corner[2]),
            new BABYLON.Vector3(corner[0], corner[1], corner[2])
        ],
        [new BABYLON.Vector3(corner[0], origin[1], corner[2]),
            new BABYLON.Vector3(corner[0], origin[1], origin[2])
        ],
    ];

   const myColors = [
        [ new BABYLON.Color4(1, 0, 0, 1),
            new BABYLON.Color4(1, 0, 0, 1)
        ],

        [ new BABYLON.Color4(0, 1, 0, 1),
            new BABYLON.Color4(0, 1, 0, 1)
        ],
        [ new BABYLON.Color4(1, 1, 0, 1),
            new BABYLON.Color4(1, 1, 0, 1)
        ],
        [ new BABYLON.Color4(1, 1, 1, 1),
            new BABYLON.Color4(1, 1, 1, 1)
        ],
        [ new BABYLON.Color4(1, 1, 1, 1),
            new BABYLON.Color4(1, 1, 1, 1)
        ],
        [ new BABYLON.Color4(1, 1, 1, 1),
            new BABYLON.Color4(1, 1, 1, 1)
        ],
        [ new BABYLON.Color4(1, 1, 1, 1),
            new BABYLON.Color4(1, 1, 1, 1)
        ],
        [ new BABYLON.Color4(1, 1, 1, 1),
            new BABYLON.Color4(1, 1, 1, 1)
        ],
        [ new BABYLON.Color4(1, 1, 1, 1),
            new BABYLON.Color4(1, 1, 1, 1)
        ],
        [ new BABYLON.Color4(1, 1, 1, 1),
            new BABYLON.Color4(1, 1, 1, 1)
        ],
        [ new BABYLON.Color4(1, 1, 1, 1),
            new BABYLON.Color4(1, 1, 1, 1)
        ],
        [ new BABYLON.Color4(1, 1, 1, 1),
            new BABYLON.Color4(1, 1, 1, 1)
        ],
    ];

    const linesystem = BABYLON.MeshBuilder.CreateLineSystem("linesystem", {lines: myLines, colors:myColors}, scene);
}
