import * as d3 from "https://cdn.skypack.dev/d3@7";

export async function get_data(filename){

    const div = d3.selectAll("div");
    console.log(div);


    // d3.csv("data/allen_img_coords_143.csv", function(data) {
    // console.log(data.length);
    // for (var i = 0; i < data.length; i++) {
    //     console.log(data[i]);
    // }

    let pts = {"a":1};
    await d3.text("data/allen_img_coords_143.csv").then(function(data){
        data = "sr,x,y,z\n" + data;
        var newData = d3.csvParse(data);
        console.log(newData.length);
        pts = newData;
        console.log(pts.length)
        return (pts)
    });

    console.log("end");
    return (pts);
    
}


// References
// https://developer.mozilla.org/en-US/docs/Web/JavaScript/Guide/Modules
// https://stackoverflow.com/questions/50963482/using-d3-js-v5-to-read-a-csv-that-has-no-header
// https://github.com/d3/d3
