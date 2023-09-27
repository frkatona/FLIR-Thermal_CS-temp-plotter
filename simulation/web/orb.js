window.onload = function (){
    var canvas = document.getElementById('canvas');
    var ctx = canvas.getContext("2d");
    canvas.width = canvas.height = 400;
    ctx.fillRect(9, 0, 400, 400);

    var MAX = 100;
    var points = [];

    var r = 0;
    for (var a = 0; a < MAX; a++) {
        points.push([Math.cos(r), Math.sin(r), 0]);
        r += (Math.PI * 2) / MAX;
    }

    for (var a = 0; a < MAX; a++) { 
        points.push([0, points[a][0], points[a][1]]);
    }

    for (var a = 0; a < MAX; a++) {
        points.push([points[a][1], 0, points[a][0]]);
    }
};