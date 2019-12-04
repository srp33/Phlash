<template>
    <!-- <div id="file-drag-drop">
        <form ref="fileform">
            <span class="drop-files">Drop the files here!</span>
        </form>
        <progress max="100" :value.prop="uploadPercentage"></progress>

        <div v-for="(file, key) in files" :key=key class="file-listing">
            <img class="preview" v-bind:ref="'preview'+parseInt( key )"/>
            {{ file.name }}
            <div class="remove-container">
                <a class="remove" v-on:click="removeFile( key )">Remove</a>
            </div>
        </div>
        <a class="submit-button" v-on:click="submitFiles()" v-show="files.length > 0">Submit</a>
    </div> -->
    <form class="box" method="post" action="" enctype="multipart/form-data">
        <div class="box__input">
            <input class="box__file" type="file" name="files[]" id="file" data-multiple-caption="{count} files selected" multiple />
            <label for="file"><strong>Choose a file</strong><span class="box__dragndrop"> or drag it here</span>.</label>
            <button class="box__button" type="submit">Upload</button>
        </div>
        <div class="box__uploading">Uploading&hellip;</div>
        <div class="box__success">Done!</div>
        <div class="box__error">Error! <span></span>.</div>
    </form>
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
    methods: {
        determineDragAndDropCapable(){
            // Create a test element to see if certain events are present that let us do drag and drop.
            var div = document.createElement('div');
            return (('draggable' in div)
                    || ('ondragstart' in div && 'ondrop' in div))
                    && 'FormData' in window
                    && 'FileReader' in window;
        },
        getImagePreviews(){
            for(let i = 0; i < this.files.length; i++) {
                if (/\.(jpe?g|png|gif)$/i.test(this.files[i].name)) {
                    let reader = new FileReader();

                    reader.addEventListener("load", function(){
                        this.$refs['preview'+parseInt( i )][0].src = reader.result;
                    }.bind(this), false);

                    reader.readAsDataURL(this.files[i]);
                } else {
                    this.$nextTick(function(){
                        this.$refs['preview'+parseInt(i)][0].src = '/images/file.png';
                    });
                }
            }
        },
        removeFile(key){
            this.files.splice( key, 1 );
        },
        submitFiles(){
            let formData = new FormData();

            for( var i = 0; i < this.files.length; i++ ){
                let file = this.files[i];

                formData.append('files[' + i + ']', file);
            }

            axios.post( '/file-drag-drop', formData, {
                headers: {
                    'Content-Type': 'multipart/form-data'
                },
                onUploadProgress: function( progressEvent ) {
                    this.uploadPercentage = parseInt( Math.round( ( progressEvent.loaded * 100 ) / progressEvent.total ) );
                }.bind(this)
            }).then(function(){
                console.log('SUCCESS!!');
            })
            .catch(function(){
                console.log('FAILURE!!');
            });
        },
    },
    mounted(){
        // Determine if drag and drop functionality is capable in the browser
        this.dragAndDropCapable = this.determineDragAndDropCapable();

        // If drag and drop capable, then we continue to bind events to our elements.
        if( this.dragAndDropCapable ){
            ['drag', 'dragstart', 'dragend', 'dragover', 'dragenter', 'dragleave', 'drop'].forEach( function( evt ) {
                this.$refs.fileform.addEventListener(evt, function(e){
                    e.preventDefault();
                    e.stopPropagation();
                }.bind(this), false);
            }.bind(this));

            this.$refs.fileform.addEventListener('drop', function(e){
                for( let i = 0; i < e.dataTransfer.files.length; i++ ){
                    this.files.push( e.dataTransfer.files[i] );
                    this.getImagePreviews();
                }
            }.bind(this));
        }
    },
}
</script>

<style>
.box__dragndrop,
.box__uploading,
.box__success,
.box__error {
  display: none;
}

.box.has-advanced-upload {
  background-color: white;
  outline: 2px dashed black;
  outline-offset: -10px;
}
.box.has-advanced-upload .box__dragndrop {
  display: inline;
}

.box.is-dragover {
  background-color: grey;
}
/*-------------------------------------*/
/* form{
    display: block;
    height: 250px;
    width: 400px;
    background: #ccc;
    margin: auto;
    margin-top: 40px;
    text-align: center;
    line-height: 400px;
    border-radius: 4px;
}

div.file-listing{
    width: 400px;
    margin: auto;
    padding: 10px;
    border-bottom: 1px solid #ddd;
}

div.file-listing img{
    height: 100px;
}

div.remove-container{
  text-align: center;
}

div.remove-container a{
  color: red;
  cursor: pointer;
}

a.submit-button{
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

progress{
  width: 400px;
  margin: auto;
  display: block;
  margin-top: 20px;
  margin-bottom: 20px;
} */
</style>