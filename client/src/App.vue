<template>
<div id="app">
   <nav class="navbar sticky-top navbar-expand-lg navbar-light bg-light">
      <a class="navbar-brand" href="/">Phlash</a>
      <button class="navbar-toggler" type="button" data-toggle="collapse" data-target="#navbarNav" aria-controls="navbarNav" aria-expanded="false" aria-label="Toggle navigation">
         <span class="navbar-toggler-icon"></span>
      </button>
      <div class="collapse navbar-collapse" id="navbarNav">
         <ul class="navbar-nav">
            <li class="nav-item">
               <router-link to="/">
                  <a class="nav-link" href="#">home</a>
               </router-link>
            </li>
            <li class="nav-item">
               <router-link :to="{ name: 'Upload', params: {currentUser: this.$route.params.currentUser} }">
                  <a class="nav-link" href="#">upload</a>
               </router-link>
            </li>
            <li class="nav-item">
               <router-link :to="{ name: 'DNAMaster', params: {currentUser: this.$route.params.currentUser} }">
                  <a class="nav-link" href="#">dna master</a>
               </router-link>
            </li>
            <li class="nav-item">
               <router-link :to="{ name: 'Blast', params: {currentUser: this.$route.params.currentUser} }">
                  <a class="nav-link" href="#">blast</a>
               </router-link>
            </li>
            <li class="nav-item">
               <router-link :to="{ name: 'Annotations', params: {currentUser: this.$route.params.currentUser} }">
                  <a class="nav-link" href="#">annotations</a>
               </router-link>
            </li>
         </ul>
      </div>
   </nav>

   <router-view/>
</div>
</template>

<style>
#app {
  font-family: 'Avenir', Helvetica, Arial, sans-serif;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
  text-align: center;
  color: #2c3e50;
  margin: 0px auto;
}
</style>
