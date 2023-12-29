class CodeWide extends HTMLElement {
    constructor() {
        super();
    }

    connectedCallback() {
        const copy = document.createElement('div');
        copy.classList.add('copy-code');
        const copy_front = document.createElement('div');
        copy_front.classList.add('front');
        const copy_back = document.createElement('div');
        copy_back.classList.add('back');
        const description = document.createElement('span');
        description.innerText = 'Copy';

        copy.addEventListener('click', async () => {
            await navigator.clipboard.writeText(this.innerText);
        });

        copy.appendChild(copy_front);
        copy.appendChild(copy_back);
        copy.appendChild(description);
        this.appendChild(copy);
    }
}

// Registra il nuovo elemento personalizzato
customElements.define('code-wide', CodeWide);
