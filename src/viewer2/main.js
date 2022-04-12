
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
    particle.position = new BABYLON.Vector3(pts[s]["x"], pts[s]["y"], pts[s]["z"]);
    // particle.position = new BABYLON.Vector3(0.5+0.25*Math.random(), i/5000, 0.25*Math.random());
    particle.color = new BABYLON.Color4(Math.random(), Math.random(), Math.random(), Math.random());
}

let makeFunc = function(inps){
    var pts = inps;
    function myFunc2(particle, i, s){
        particle.position = new BABYLON.Vector3(pts[s]["x"], pts[s]["y"], pts[s]["z"]);
        // particle.position = new BABYLON.Vector3(0.5+0.25*Math.random(), i/5000, 0.25*Math.random());
        particle.color = new BABYLON.Color4(Math.random(), Math.random(), Math.random(), Math.random());

    }
    return myFunc2;

}

let sliceSelector = function(state){

    let showState = state;
    function updateParticle(particle){
        if (showState==true){
            particle.color = new BABYLON.Color4(1.0, 0, 0, 1.0);
        }else{
            particle.color = new BABYLON.Color4(0.8, 0.8, 0.8, 0.0);

        }
    }
    return updateParticle;
}

const zeroPad = (num, places) => String(num).padStart(places, '0');

export async function create_babylon () {

    let data_map = new Map();

    const nis_ids = [];
    const avoid_list = [5, 77, 167, 181, 205, 223, 225, 227];
    for (let i=1; i<50; i++){
        if (!avoid_list.includes(i) && (i%2==1))
            nis_ids.push(i);
        
    }
    console.log(nis_ids);

    const get_data_promises = [];
    // for (let i=141; i<144; i+=2){
    // for (const i of nis_ids){
    nis_ids.forEach(function (nisid, idx) {
        let zpadded_i = zeroPad(nisid, 3);
        let file = `data/allen_img_coords_${zpadded_i}.csv`;
        get_data_promises.push(get_data(file));
    });

    let all_slice_data = await Promise.all(get_data_promises);

    nis_ids.forEach(function (nisid, idx) {
        console.log('%d: %s', idx, nisid);
        data_map.set(nisid, all_slice_data[idx] )
    });

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

        const camera = new BABYLON.ArcRotateCamera("camera", -Math.PI / 2, Math.PI / 2.5, 1200, new BABYLON.Vector3(250, 250, 250));
        camera.attachControl(canvas, true);

        const light = new BABYLON.HemisphericLight("light", new BABYLON.Vector3(0, 1, 0));

        // const box = BABYLON.MeshBuilder.CreateBox("box", {});

        let pcs_map = new Map();
        let point_size = 2;
        var pcs = new BABYLON.PointsCloudSystem("pcs", point_size, scene) 
        nis_ids.forEach((nisid, idx) => {
            pcs_map.set(nisid, new BABYLON.PointsCloudSystem(`pcs{nisid}`, point_size, scene));  
            let pts = data_map.get(nisid);
            console.log(pts);
            let myFunc = makeFunc(pts);
            let pcs = pcs_map.get(nisid);
            // console.log(pcs.mesh.material);
            pcs.addPoints(pts.length, myFunc);
            pcs.buildMeshAsync().then(()=>{
                console.log(pcs.mesh.visibility);
                // pcs.mesh.visibility = 0.1;

            });
        });

        // pcs.updateParticle = function(particle){
        //     particle.rotation.y +=0.1;
        // }

        let selected_slice = -1; // modify here to jump directly to slice
        let next_slice = -1;
        let nSlices = nis_ids.length;

        scene.onKeyboardObservable.add((kbInfo) => {
            if (kbInfo.event.key==='j')
                next_slice = (selected_slice+2)%(nSlices*2);
            else if (kbInfo.event.key==='k'){
                next_slice = (selected_slice-2)%(nSlices*2);
                if (next_slice<0)
                    next_slice = (nSlices*2-1);
            }
            switch (kbInfo.type) {
                case BABYLON.KeyboardEventTypes.KEYDOWN:
                    // console.log("KEY DOWN: ", kbInfo.event.key);
                    selected_slice = next_slice;
                    nis_ids.forEach((nisid, idx) => {
                        pcs = pcs_map.get(nisid);
                        if (nisid==selected_slice){
                            pcs.updateParticle = sliceSelector(true);
                            pcs.mesh.visibility = 1.0;
                        }else{
                            pcs.updateParticle = sliceSelector(false);
                            pcs.mesh.visibility = 0.02;
                        }
                        pcs.setParticles();
                    });
                    console.log("KEY DOWN: ", kbInfo.event.key==='j', "selected_slice", selected_slice, "nSlices*2", nSlices*2, "selected_slice", selected_slice);
                    break;
                // case BABYLON.KeyboardEventTypes.KEYUP:
                //     console.log("KEY UP: ", kbInfo.event.code);
                //     break;
            }
        });

        return scene;
    }

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
// https://stackoverflow.com/questions/2998784/how-to-output-numbers-with-leading-zeros-in-javascript
