<template>
   <div id="file-drag-drop">
      <form ref="fileform" type="file" enctype="multipart/form-data">
         <span class="drop-files">Drop the files here or</span><br>
         <label id="label-file" for="file">Select files:</label>
         <input class="browse-files" type="file" name="files[]" id="file" data-multiple-caption="{count} files selected" multiple /><br>
      </form>

      <progress max="100" :value.prop="uploadPercentage"></progress>

      <div v-for="(file, key) in files" :key="key" class="file-listing">
         <img class="preview" v-bind:ref="'preview'+parseInt( key )"/>
         {{ file.name }}
         <div class="remove-container">
            <a class="remove" v-on:click="removeFile( key )">Remove</a>
         </div>
      </div>

      <button class="btn btn-primary" @click="submitFiles" v-show="files.length > 0">Upload</button>
   </div>
</template>

<script>
import axios from 'axios';

export default {
   data() {
      return {
         dragAndDropCapable: false,
         files: [],
         uploadPercentage: 0
      }
   },
   mounted(){
      this.dragAndDropCapable = this.determineDragAndDropCapable();

      if( this.dragAndDropCapable ){

         ['drag', 'dragstart', 'dragend', 'dragover', 'dragenter', 'dragleave', 'drop'].forEach( function( evt ) {
            this.$refs.fileform.addEventListener(evt, function(e) {
               e.preventDefault();
               e.stopPropagation();
            }.bind(this), false);
         }.bind(this));

         this.$refs.fileform.addEventListener('drop', function(e) {
            // this.files = this.$refs.fileform.files

            for( let i = 0; i < e.dataTransfer.files.length; i++ ) {
               this.files.push( e.dataTransfer.files[i] );
            }
         }.bind(this));
      }
   },
   methods: {
      determineDragAndDropCapable() {
         var div = document.createElement('div');
         return ( ( 'draggable' in div )
                  || ( 'ondragstart' in div && 'ondrop' in div ) )
                  && 'FormData' in window
                  && 'FileReader' in window;
      },
      submitFiles(){
         let formData = new FormData();
         for( var i = 0; i < this.files.length; i++ ) {
            let file = this.files[i];

            formData.append('files[' + i + ']', file);
         }
         axios.post('http://localhost:5000/api/upload',
            formData,
            {
               headers: {
                  'Content-Type': 'multipart/form-data'
               },
               onUploadProgress: function( progressEvent ) {
                  this.uploadPercentage = parseInt( Math.round( ( progressEvent.loaded * 100 ) / progressEvent.total ) );
               }.bind(this)
            }
         ).then(function(){
            console.log('SUCCESS!!');
         })
         .catch(function(){
            console.log('FAILURE!!');
         });
      },
      removeFile( key ){
         this.files.splice( key, 1 );
      },
   }
}
</script>

<style scoped>
/* .box__dragndrop,
.box__uploading,
.box__success,
.box__error {
  display: none;
} */

/*--------------------------------*/
form {
   display: block;
   height: 200px;
   width: 400px;
   background: #ccc;
   margin: auto;
   margin-top: 40px;
   /* text-align: center; */
   /* line-height: 400px; */
   border-radius: 4px;
   position: relative;
}

.label-file {
   position: absolute;
   transform: translate(-50%, 150%);
}

.drop-files {
   position: absolute;
   transform: translate(-50%, 130%);
}

.browse-files {
   position: absolute;
   transform: translate(-40%, 200%);
}

div.file-listing {
   width: 400px;
   margin: auto;
   padding: 10px;
   border-bottom: 1px solid #ddd;
}

div.file-listing img {
   height: 100px;
}

div.remove-container {
   text-align: center;
}

div.remove-container a {
   color: red;
   cursor: pointer;
}

a.submit-button {
   display: block;
   margin: auto;
   text-align: center;
   width: 200px;
   padding: 10px;
   text-transform: uppercase;
   background-color: #CCC;
   color: white;
   font-weight: bold;
   margin-top: 20px;
}

progress {
   width: 400px;
   margin: auto;
   display: block;
   margin-top: 20px;
   margin-bottom: 20px;
}
</style>