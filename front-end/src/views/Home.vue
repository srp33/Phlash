<template>
<div class="container">
   <h1><strong>Phlash</strong></h1>

   <div class="alert alert-primary" align="left">
      <p>
         Welcome to <strong><i>Phlash</i></strong>! Enter an ID below to get started.
      </p>
      <div class="input-group mb-3">
         <input type="text" v-model="id" class="form-control" v-on:keyup.enter="checkID(id)" placeholder="Enter a unique ID" aria-label="Enter a unique ID" aria-describedby="basic-addon2">
         <div class="input-group-append">
            <button class="btn btn-dark btn-sm" type="button" @click="checkID(id)">
               <strong>Enter</strong>
            </button>
         </div>
      </div>
      <p class="message">{{ message }}</p>
      <router-link :to="{ name: 'DNAMaster', params: {currentUser: id} }" v-show="message.includes('ID already exists')">
         <button class="btn btn-light"><strong>Next</strong></button>
      </router-link>
      <router-link :to="{ name: 'Upload', params: {currentUser: id} }" v-show="message.includes('ID created')">
         <button class="btn btn-light"><strong>Next</strong></button>
      </router-link>
   </div>
</div>
</template>

<script>
import axios from 'axios';

export default {
   name: 'Home',
   data() {
      return {
         id: "",
         message: "",
      }
   },
   methods: {
      checkID(id) {
         axios.post(`http://localhost:5000/api/home/${id}`)
         .then(response => {
            this.message = response.data.message
            // this.$root.$data.currentUser = this.id
            console.log(this.id)
         })
         .catch(error => {
            console.error(error)
         });
      }
   }
}
</script>

<style scoped>
h1 {
   margin: 40px auto;
}

.message {
   margin-top: 10px;
   font-style: italic;
}
</style>
