<template>
    <div class="container">
        <h3>Upload FASTA and DNA Master Files</h3>
        <br>
        <div align="left">
            <!-- <h2>Docs</h2>
            <ul>
                <li>What does this app do?</li>
                <li>DNAMaster instructions</li>
                <li>GeneMark instructions</li>
            </ul>
            <p>fixme: add alert for successfully uploading file</p> -->
            <br>
            <p>Please upload a FASTA file and your genbank file outputted from DNA Master.</p>
            <br>
        </div>
        <alert :message="message" v-if="showMessage"></alert>
        <form id="upload_form" role="form" enctype="multipart/form-data">
            <input type="file" id="files" ref="files" multiple v-on:change="handleFilesUpload" class="form-control">
            <button class="btn btn-success btn-block" @click="upload">Upload</button>
        </form>
        <br>
        <button class="btn btn-primary" @click="updateDatabase">Add DNA Master data to database</button>
        <br><br>
        <router-link  to="/database" v-if="showDatabase">
                <button class="btn btn-success"><strong>Next: View your data in the database</strong></button>
        </router-link>
    </div>
</template>

<script>
import axios from 'axios';
import Alert from './Alert.vue';

export default {
    data() {
        return {
            files: null,
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
            axios.post('http://localhost:5000/api/upload',
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
            axios.get('http://localhost:5000/api/upload')
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
