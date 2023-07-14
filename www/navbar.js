nav = document.querySelector(".navbar-header")
navContents = document.querySelector("nav .container-fluid")
navContents.classList.add("collapsed-nav")

burger = document.createElement('input');
burger.type = 'checkbox';
burger.id = "burger"
burger.checked = false;

burgerLabel = document.createElement("label")
burgerLabel.htmlFor = "burger"
burgerLabel.id = "burger-label"

nav.appendChild(burger);
nav.appendChild(burgerLabel);

burger.addEventListener("change", (e) => {
    console.log(burger.checked)
    burger.checked ? navContents.classList.remove("collapsed-nav") : navContents.classList.add("collapsed-nav")     
})
