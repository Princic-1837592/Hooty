const eyes = document.getElementsByClassName("hooty-eye-ball") as HTMLCollectionOf<HTMLElement>;
const eyeCenter = document.getElementById("hooty-eye-center");

document.addEventListener("mousemove", (e) => {
    const mouseX = e.clientX;
    const mouseY = e.clientY;
    const eyeCenterRect = eyeCenter.getBoundingClientRect();
    const eyeCenterX = eyeCenterRect.left + eyeCenterRect.width / 2;
    const eyeCenterY = eyeCenterRect.top + eyeCenterRect.height / 2;
    for (let i = 0; i < eyes.length; i++) {
        const eye = eyes[i];
        const [dx, dy, d] = distance(eyeCenterX, eyeCenterY, mouseX, mouseY);
        if (d < 50) {
            eye.classList.add("staring");
            eye.classList.remove("rotating");
        } else {
            const deg = angle(dx, dy);
            eye.style.transform = `rotate(${deg}deg)`;
            eye.classList.add("rotating");
            eye.classList.remove("staring");
        }
    }
});

function distance(cx: number, cy: number, ex: number, ey: number): [number, number, number] {
    const dx = ex - cx;
    const dy = ey - cy;
    return [dx, dy, Math.sqrt(Math.pow(dx, 2) + Math.pow(dy, 2))];
}

function angle(dx: number, dy: number): number {
    const rad = Math.atan2(dy, dx);
    return rad * 180 / Math.PI;
}
