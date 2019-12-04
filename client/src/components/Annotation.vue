<template>
    <div class="container">
        <h1>Annotations</h1>
        <br>
        <div id="annotations" align="center">
            <div align="left">
                <h5>Status</h5>
                <p>
                    <strong style="color:green">Pass:</strong> DNA Master's gene call is the same as or longer than GeneMark's gene call.<br />
                    <strong style="color:red">Fail:</strong> DNA Master's gene call is shorter than GeneMark's gene call.<br /> 
                    <strong style="color:orange">Not called by GeneMark:</strong> DNA Master's gene call has not been called at all in GeneMark.<br/>
                </p>
                <p>
                    The "Action" column must contain only "Done"s for you to continue. Once all anotations are complete, click the button below to continue.
                </p>
            </div>
            <button class="btn btn-primary" disabled><strong>Create GenBank File</strong></button>
            <br><br>
            <table class="table table-hover" align="center">
                <thead>
                    <tr>
                        <th scope="col">ID</th>
                        <th scope="col">Start</th>
                        <th scope="col">Stop</th>
                        <th scope="col">Strand</th>
                        <th scope="col">Function</th>
                        <th scope="col">Status</th>
                        <th scope="col">Action</th>
                    </tr>
                </thead>
                <tbody>
                    <tr v-for="(curr, index) in dnamaster" :key="index">
                        <td>{{ curr.id }}</td>
                        <td>{{ curr.start }}</td>
                        <td>{{ curr.stop }}</td>
                        <td>{{ curr.strand }}</td>
                        <td v-if="curr.function.length<30">{{ curr.function }}</td>
                        <td v-else>{{ curr.function.substring(0,30) }}...</td>
                        <td>{{ curr.status }}</td>
                        <td>
                            <span style="color:green" 
                                  v-if="curr.function !== 'None' && curr.status == 'Pass'">
                                <strong>Done</strong>
                            </span>
                            <router-link :to="{ name: 'Blast', params: {id: curr.id} }">
                                <button class="btn btn-success btn-sm" 
                                        v-if="curr.function == 'None' && curr.status == 'Pass'">
                                    <strong>BLAST for function</strong>
                                </button>
                            </router-link>
                            <router-link :to="{ name: 'Failed', params: {id: curr.id} }">
                                <button class="btn btn-danger btn-sm" style="width:141px"
                                        v-if="curr.status == 'Fail'">
                                    <strong>Modify gene call</strong>
                                </button>
                            </router-link>
                            <router-link :to="{ name: 'More', params: {id: curr.id} }">
                                <button class="btn btn-warning btn-sm" style="width:141px"
                                        v-if="curr.status =='Not called by GeneMark'">
                                    <strong>See more</strong>
                                </button>
                            </router-link>
                        </td>
                    </tr>
                </tbody>
            </table>
        </div>
    </div>
</template>

<script>
import axios from 'axios';

export default {
    data() {
        return {
            dnamaster: [],
            showFail: false,
            updatedCDS: {
                id: '',
                start: '',
                stop: '',
                strand: '',
                read: [],
            },
        };
    },
    methods: {
        getData() {
            axios.get('http://localhost:5000/annotate_data')
            .then(response => {
                this.dnamaster = response.data.dnamaster;
            })
            .catch(error => {
                console.error(error);
            });
        },
    },
    created() {
        this.getData();
    },
};
</script>

<style scoped>
.button {
    width: 250px;
    text-align: center;
}
</style>