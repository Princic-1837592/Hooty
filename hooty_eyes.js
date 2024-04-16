var eyes = document.getElementsByClassName("hooty-eye-ball");
var eyeCenter = document.getElementById("hooty-eye-center");
document.addEventListener("mousemove", function (e) {
    var mouseX = e.clientX;
    var mouseY = e.clientY;
    var eyeCenterRect = eyeCenter.getBoundingClientRect();
    var eyeCenterX = eyeCenterRect.left + eyeCenterRect.width / 2;
    var eyeCenterY = eyeCenterRect.top + eyeCenterRect.height / 2;
    for (var i = 0; i < eyes.length; i++) {
        var eye = eyes[i];
        var _a = distance(eyeCenterX, eyeCenterY, mouseX, mouseY), dx = _a[0], dy = _a[1], d = _a[2];
        if (d < 50) {
            eye.classList.add("staring");
            eye.classList.remove("rotating");
        }
        else {
            var deg = angle(dx, dy);
            eye.style.transform = "rotate(".concat(deg, "deg)");
            eye.classList.add("rotating");
            eye.classList.remove("staring");
        }
    }
});
function distance(cx, cy, ex, ey) {
    var dx = ex - cx;
    var dy = ey - cy;
    return [dx, dy, Math.sqrt(Math.pow(dx, 2) + Math.pow(dy, 2))];
}
function angle(dx, dy) {
    var rad = Math.atan2(dy, dx);
    return rad * 180 / Math.PI;
}
