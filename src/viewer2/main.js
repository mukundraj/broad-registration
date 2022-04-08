
import {get_data} from './loaders.js'
// import 'https://cdnjs.cloudflare.com/ajax/libs/dat-gui/0.6.2/dat.gui.min.js'
// import "https://assets.babylonjs.com/generated/Assets.js"
import "https://preview.babylonjs.com/ammo.js"
// import "https://preview.babylonjs.com/cannon.js"
// import "https://preview.babylonjs.com/Oimo.js"
// import "https://preview.babylonjs.com/earcut.min.js"
import "https://preview.babylonjs.com/babylon.js"
// import "https://preview.babylonjs.com/materialsLibrary/babylonjs.materials.min.js"
// import "https://preview.babylonjs.com/proceduralTexturesLibrary/babylonjs.proceduralTextures.min.js"
// import "https://preview.babylonjs.com/postProcessesLibrary/babylonjs.postProcess.min.js"
// import "https://preview.babylonjs.com/loaders/babylonjs.loaders.js"
// import "https://preview.babylonjs.com/serializers/babylonjs.serializers.min.js"
// import "https://preview.babylonjs.com/gui/babylon.gui.min.js"
// import "https://preview.babylonjs.com/inspector/babylon.inspector.bundle.js"

export function main(){

    console.log("in main");
    get_data("data/allen_img_coords_143.csv")

    console.log("done")


    var canvas = document.getElementById("renderCanvas");

    var startRenderLoop = function (engine, canvas) {
        engine.runRenderLoop(function () {
            if (sceneToRender && sceneToRender.activeCamera) {
                sceneToRender.render();
            }
        });
    }

    var engine = null;
    var scene = null;
    var sceneToRender = null;
    var createDefaultEngine = function() { return new BABYLON.Engine(canvas, true, { preserveDrawingBuffer: true, stencil: true,  disableWebGL2Support: false}); };

    const createScene =  () => {
        const scene = new BABYLON.Scene(engine);

        const camera = new BABYLON.ArcRotateCamera("camera", -Math.PI / 2, Math.PI / 2.5, 3, new BABYLON.Vector3(0, 0, 0));
        camera.attachControl(canvas, true);

        const light = new BABYLON.HemisphericLight("light", new BABYLON.Vector3(0, 1, 0));

        const box = BABYLON.MeshBuilder.CreateBox("box", {});

        var pcs= new BABYLON.PointsCloudSystem("pcs", 5, scene) 
        //pcs.mesh.material.pointSize = 1;

        pcs.addPoints(10000);
        pcs.buildMeshAsync();

        return scene;
    }
    window.initFunction = async function() {


        var asyncEngineCreation = async function() {
            try {
                return createDefaultEngine();
            } catch(e) {
                console.log("the available createEngine function failed. Creating the default engine instead");
                return createDefaultEngine();
            }
        }

        window.engine = await asyncEngineCreation();
        if (!engine) throw 'engine should not be null.';
        startRenderLoop(engine, canvas);
        window.scene = createScene();};
    initFunction().then(() => {sceneToRender = scene                    
    });

    // Resize
    window.addEventListener("resize", function () {
        engine.resize();
    });
}


let myFunc = function (particle, i, s){
    particle.position = new BABYLON.Vector3(pts[i]["x"], pts[i]["y"], pts[i]["z"]);
    // particle.position = new BABYLON.Vector3(0.5+0.25*Math.random(), i/5000, 0.25*Math.random());
    particle.color = new BABYLON.Color4(Math.random(), Math.random(), Math.random(), Math.random());
}

let makeFunc = function(inps){
    var pts = inps;
    function myFunc2(particle, i, s){
        particle.position = new BABYLON.Vector3(pts[i]["x"], pts[i]["y"], pts[i]["z"]);
        // particle.position = new BABYLON.Vector3(0.5+0.25*Math.random(), i/5000, 0.25*Math.random());
        particle.color = new BABYLON.Color4(Math.random(), Math.random(), Math.random(), Math.random());

    }
    return myFunc2;

}


export async function create_babylon () {
    let pts = await get_data("data/allen_img_coords_143.csv")
    console.log(pts)
    var canvas = document.getElementById('renderCanvas')

    var sceneToRender = null
    var createDefaultEngine = function () {
        return new BABYLON.Engine(canvas, true, {
            preserveDrawingBuffer: true,
            stencil: true,
            disableWebGL2Support: false
        })
    }

    var createDefaultEngine = function(engine) { return new BABYLON.Engine(canvas, true, { preserveDrawingBuffer: true, stencil: true,  disableWebGL2Support: false}); };

    const createScene =  () => {
        const scene = new BABYLON.Scene(engine);

        const camera = new BABYLON.ArcRotateCamera("camera", -Math.PI / 2, Math.PI / 2.5, 3, new BABYLON.Vector3(0, 0, 0));
        camera.attachControl(canvas, true);

        const light = new BABYLON.HemisphericLight("light", new BABYLON.Vector3(0, 1, 0));

        // const box = BABYLON.MeshBuilder.CreateBox("box", {});

        var pcs = new BABYLON.PointsCloudSystem("pcs", 2, scene) 
        //pcs.mesh.material.pointSize = 1;

        // pcs.addPoints(10000);
        let myFunc3 = makeFunc(pts);
        pcs.addPoints(pts.length, myFunc3)
        pcs.buildMeshAsync();

        return scene;
    }
    // const createScene = (engine) => {
    //     const scene = new BABYLON.Scene(engine)

    //     const camera = new BABYLON.ArcRotateCamera('Camera', 3 * Math.PI / 4, Math.PI / 4, 4, BABYLON.Vector3.Zero(), scene)
    //     camera.attachControl()
    //     const light = new BABYLON.HemisphericLight(
    //         'light',
    //         new BABYLON.Vector3(1, 1, 0),
    //         scene, // Always pass this argument explicitly
    //     )

    //     const box = BABYLON.MeshBuilder.CreateBox(
    //         'box',
    //         {height: 1, width: 0.75, depth: 0.25},
    //         scene, // Always pass this argument explicitly
    //     )

    //     return scene

    // }

    window.initFunction = async function () {

        var asyncEngineCreation = function () {
            try {
                return createDefaultEngine()
            } catch (e) {
                console.log('the available createEngine function failed. Creating the default engine instead')
                return createDefaultEngine()
            }
        }

        window.engine = asyncEngineCreation()
        if (!window.engine) throw 'engine should not be null.'
        window.scene = createScene(window.engine)
    }
    initFunction().then(() => {
        sceneToRender = scene
        engine.runRenderLoop(function () {
            if (sceneToRender && sceneToRender.activeCamera) {
                sceneToRender.render()
            }
        })
    })

    // Resize
    window.addEventListener('resize', function () {
        engine.resize()
    })

}

// References
// https://forum.babylonjs.com/t/create-babylon-instance-with-a-function-uncaught-in-promise-engine-should-not-be-null/26588/11
