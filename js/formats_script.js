var ul = document.getElementById('myUL')
var myNodelist = ul.getElementsByTagName("LI");
var i;
var smiles_file_data="";
var txt = document.createTextNode("\u00D7");
for (i = 0; i < myNodelist.length; i++) {
  var span = document.createElement("SPAN");
  var txt = document.createTextNode("\u00D7");
  span.className = "close";
  span.appendChild(txt);
  myNodelist[i].appendChild(span);
}

var close = document.getElementsByClassName("close");
var i;
for (i = 0; i < close.length; i++) {
  close[i].onclick = function() {
    var div = this.parentElement;
    div.style.display = "none";
  }
}

// Add a "checked" symbol when clicking on a list item
var list = document.querySelector('.listSmiles');
list.addEventListener('click', function(ev) {
  if (ev.target.tagName === 'LI') {
    document.querySelector('.checked').classList.remove('checked');
    ev.target.classList.toggle('checked');
    let li = document.querySelector('[class*="checked"]');
    let smile= li.textContent;
    showFormat(smile);
    // resetViewer()
    // get3D(smile)
  }
}, false);

// Create a new list item when clicking on the "Add" button
function newElement() {
  // resetViewer()
  // readSmileFromInput();
  var li = document.createElement("li");
  var inputValue = document.getElementById("myInput").value;
  var t = document.createTextNode(inputValue);
  showFormat(inputValue);
  li.appendChild(t);
  if (inputValue === '') {
    alert("You must write something!");
  } else {
    document.getElementById("myUL").appendChild(li);
  }
  document.getElementById("myInput").value = "";

  li.setAttribute('id','smileLi');
  document.querySelector('.checked').classList.remove('checked');
  li.classList.toggle('checked');

  var span = document.createElement("SPAN");
  var txt = document.createTextNode("\u00D7");
  span.className = "close";
  span.appendChild(txt);
  li.appendChild(span);

  for (i = 0; i < close.length; i++) {
    close[i].onclick = function() {
      var div = this.parentElement;
      div.style.display = "none";
    }
  }

}

document.getElementById('fileToLoad').onchange = function(){

  var file = this.files[0];

  var reader = new FileReader();
  reader.onload = function(progressEvent){
    console.log(this.result);
    if (window.File && window.FileReader && window.FileList && window.Blob) {
        var file = document.querySelector('input[type=file]').files[0];
        var reader = new FileReader()

        reader.onload = function (event) {
              $("#textBox-smile").val(event.target.result);
        }
        reader.readAsText(file);
    }
    var lines = this.result.split('\n');
    // smiles_file_data = lines;
    for(var line = 0; line < lines.length; line++){
      console.log("Smile:"+lines[line]);
      var li = document.createElement("li");
      var inputValue = lines[line];
      var t = document.createTextNode(inputValue);
      li.appendChild(t);
      if (inputValue === '') {
        alert("You must write something!");
      } else {
        document.getElementById("myUL").appendChild(li);
      }
      document.getElementById("myInput").value = "";

      var span = document.createElement("SPAN");
      var txt = document.createTextNode("\u00D7");
      span.className = "close";
      span.appendChild(txt);
      li.appendChild(span);


      for (i = 0; i < close.length; i++) {
        close[i].onclick = function() {
          var div = this.parentElement;
          div.style.display = "none";
        }
      }
    }

  };
  reader.readAsText(file);
};

function get3D(smile){
  $.get('http://127.0.0.1:8000/3d/'+smile, function(data) {
  let cell_viewer = window.viewers[0][0];
  cell_viewer.addModel(data,'sdf');
  cell_viewer.setStyle({stick:{}});
  // cell_viewer.addSurface($3Dmol.SurfaceType.VDW, {
  //               opacity:0.85,
  //               voldata: new $3Dmol.VolumeData(volumedata, "cube"),
  //               volscheme: new $3Dmol.Gradient.ROYGB(range[1],range[0])
  //           },{});
  cell_viewer.zoomTo();     
  cell_viewer.render();

  // cell_viewer = window.viewers[0][1];
  // cell_viewer.addModel(data,'sdf');
  // cell_viewer.setStyle({sphere:{}});
  // cell_viewer.zoomTo();
  // cell_viewer.render();


});
}
function resetViewer(){
  window.viewers = $3Dmol.createViewerGrid(
    '3dmol-div', //id of div to create canvas in
    {
      rows:1,
      cols:1,
      control_all: true  //mouse controls all viewers
    },
    { backgroundColor: 'black' }
  );
}

function getFormatDoc(smile,format){
  let responseData = null;
  $.get({
    url: 'http://127.0.0.1:8000/format/'+smile+"_"+format,
    async:false,
    success: function(response) {
      console.log('http://127.0.0.1:8000/format/', response);
      responseData = response.data;
    }
  });
  return responseData;
}

function getSmiles_Zip(){
  let ul = document.getElementById('myUL')
  let listItems = ul.getElementsByTagName("LI");
  // let span = document.querySelector("close");
  let i;
  let smiles_file="";
  for (i = 0; i < listItems.length; i++) {
      // let elem = listItems[i].removeChild(listItems[i].lastChild);
      smiles_file+=listItems[i].textContent.replace('\u00D7',"")+",";

  }
  console.log(smiles_file);
  return smiles_file;

}

function showFormat(smile){
  let format =$("#formatsDropdown").val();
  let textFormat = getFormatDoc(smile,format);
  console.log("Showing showFormat()", textFormat);
  $("#textBox").val(textFormat);
}

function downloadFormat(){
  let li = document.querySelector('[class*="checked"]');
  let smile= li.textContent;
  let format =$("#formatsDropdown").val();
  let url = 'http://127.0.0.1:8000/download/'+smile+'/'+format;
  window.open(url, '_blank');
}
function getFileForZIP(){
  let format =$("#formatsDropdown").val();
  smiles=getSmiles_Zip();
  let url = 'http://127.0.0.1:8000/zip/'+smiles+'/'+format;
  window.open(url, '_blank');
}

function readSmileFromInput(){
  let smileTxt = $("input#myInput").val();
  console.log("myInput: ", smileTxt)
  if(smileTxt != "" && smileTxt != null){
    showFormat(smile);
  }else{
    alert("Input text empty.")
  }
}

document.getElementById('formatsDropdown').onchange = function () {
  let li = document.querySelector('[class*="checked"]');
  let smile= li.textContent;
  showFormat(smile);
};

$( document ).ready(function() {
  showFormat("BrC(Br)(Br)C1=NC(=NC(=N1)C(Br)(Br)Br)C(Br)(Br)Br")
});

$( ".button-group > div" ).click(function() {
  $('.button-group > div.active').not(this).removeClass('active');
  $( this ).toggleClass( "active" );
});

$('.value').each(function() {
	var text = $(this).text();
	$(this).parent().css('width', text);
});

$('.block').tooltip();

