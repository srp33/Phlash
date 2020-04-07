<template>
<div class="container">
    <h2><strong>Create GenBank File</strong></h2>
    <br>
    <div>
        <h5 align="left">Submitting author(s)</h5><br>
        <table class="table">
            <thead>
                <tr>
                    <td><strong>First name</strong></td>
                    <td><strong>M.I.</strong></td>
                    <td><strong>Last name</strong></td>
                    <td><strong>Suffix name</strong></td>
                    <td></td>
                </tr>
            </thead>
            <tbody>
                <tr v-for="(author, index) in authors" :key="index">
                    <td><input type="text" v-model="author.firstName"></td>
                    <td><input type="text" v-model="author.middleInitial"></td>
                    <td><input type="text" v-model="author.lastName"></td>
                    <td><input type="text" v-model="author.suffix"></td>
                    <td>
                        <a v-on:click="removeElement(index);" style="cursor: pointer">Remove</a>
                    </td>
                </tr>
            </tbody>
        </table>
        <div>
            <button class="button btn-primary" @click="addAuthor">Add another author</button>
        </div>
        <!-- ------------------------------------ -->
        <b-container fluid>
            <b-row class="my-2">
                <b-col sm="2" align="right">
                    <label>First name</label>
                </b-col>
                <b-col sm="8">
                    <b-form-input :id="`type-text`" :type="text"></b-form-input>
                </b-col>
            </b-row>
            <b-row class="my-2">
                <b-col sm="2" align="right">
                    <label>M.I.</label>
                </b-col>
                <b-col sm="8">
                    <b-form-input :id="`type-text`" :type="text"></b-form-input>
                </b-col>
            </b-row>
            <b-row class="my-2">
                <b-col sm="2" align="right">
                    <label>Last name</label>
                </b-col>
                <b-col sm="8">
                    <b-form-input :id="`type-text`" :type="text"></b-form-input>
                </b-col>
            </b-row>
            <b-row class="my-2">
                <b-col sm="2" align="right">
                    <label>Suffix</label>
                </b-col>
                <b-col sm="8">
                    <b-form-input :id="`type-text`" :type="text"></b-form-input>
                </b-col>
            </b-row><br>
        </b-container>
        <button align="center" class="btn btn-outline-primary btn-sm" @click="showAuthorForm">Add another author</button>
        <author v-if="showAuthor"></author>
    </div>
</div>
</template>

<script>
import axios from 'axios';

export default {
    data() {
        return {
            authors: [],
        }
    },
    components: {
    },
    methods: {
        addAuthor() {
            var elem = document.createElement('tr');
            this.authors.push({
                firstName: "",
                middleInitial: "",
                lastName: "",
                suffix: "",
            });
        },
        removeAuthort(index) {
            this.authors.splice(index, 1);
        },
    },
}
</script>