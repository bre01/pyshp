const fs = require('fs');
let file=fs.readFileSync('a.json');
console.log(JSON.parse(file));
