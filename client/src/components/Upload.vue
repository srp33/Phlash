<template>
    <div class="container">
        <h1><strong><i>PHLASH</i></strong></h1>
        <br>
        <div align="left">
            <h2>Docs</h2>
            <ul>
                <li>What does this app do?</li>
                <li>DNAMaster instructions</li>
                <li>GeneMark instructions</li>
            </ul>
            <p>fixme: add alert for successfully uploading file</p>
            <br>
            <h5>Upload your genbank file outputted from DNA Master to get started!</h5>
            <br>
        </div>
        <alert :message="message" v-if="showMessage"></alert>
        <form id="upload_form" role="form" enctype="multipart/form-data">
            <input type="file" id="file" ref="file" v-on:change="handleFileUpload" class="form-control">
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
            file: null,
            message: '',
            showMessage: false,
            showDatabase: false,
        }
    },
    components: {
        alert: Alert,
    },
    methods: {
        handleFileUpload() {
            this.file = this.$refs.file.files[0]
        },
        upload() {
            var data = new FormData();
            data.append('file', this.file);
            axios.post('http://localhost:5000/api/upload', data)
            .then(() => {
                this.message = 'File successfully uploaded!';
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
