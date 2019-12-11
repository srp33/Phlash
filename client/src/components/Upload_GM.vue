<template>
    <div class="container">
        <h3>Upload GeneMark Files</h3>
        <br>
        <div align="left">
            <!-- <p>FIXME: Try drag-and-drop multiple files uploader</p> -->
            <!-- <p>Improve docs</p> -->
            <p>Please upload the following files from GeneMark:</p>
            <ul>
                <!-- <li>.fasta.gdata (GeneMark)</li> -->
                <li>.fasta.ldata</li>
                <li>.fasta.gdata</li>
            </ul>
        </div>
        <alert :message="message" v-if="showMessage"></alert>
        <form id="upload_form" role="form" enctype="multipart/form-data">
            <input type="file" id="files" ref="files" multiple v-on:change="handleFilesUpload" class="form-control">
            <button class="btn btn-success btn-block" @click="upload">Upload</button>
        </form>
        <br>
        <button class="btn btn-primary" @click="updateDatabase">Add GeneMark data to database</button>
        <br><br>
        <router-link  to="/annotate_data" v-if="showDatabase">
                <button class="btn btn-success"><strong>Next: Do analysis!</strong></button>
        </router-link>
    </div>
</template>

<script>
import axios from 'axios';
import Alert from './Alert.vue';

export default {
    data() {
        return {
            files: '',
            message: '',
            showMessage: false,
            showDatabase: false,
        }
    },
    components: {
        alert: Alert,
    },
    methods: {
        handleFilesUpload() {
            this.files = this.$refs.files.files
        },
        upload() {
            var data = new FormData();
            // data.append('file', this.file);
            for (var i = 0; i < this.files.length; i++) {
                let file = this.files[i];
                data.append('files[' + i + ']', file);
            }
            axios.post('http://localhost:5000/api/upload_genemark',
                data,
                {
                    headers: {
                        'Content-Type': 'multipart/form-data'
                    }
                }
            ).then(() => {
                this.message = 'Files successfully uploaded!';
                this.showMessage = true;
                console.log(this.message)
            })
            .catch(error => {
                console.log(error)
            });
        },
        updateDatabase() {
            this.message = 'Updating the database will take a few seconds...';
            this.showMessage = true;
            axios.get('http://localhost:5000/api/upload_genemark')
            .then(() => {
                this.message = "Finished updating database.";
                this.showMessage = true;
                this.showDatabase = true;
            })
            .catch(error => {
            console.log(error)
            });
        },
    },
};
</script>
