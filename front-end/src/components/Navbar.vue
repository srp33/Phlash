<template>
  <div class="navbar-wrapper">
    <nav class="navbar sticky-top navbar-expand-lg navbar-light bg-light">
      <a class="navbar-brand" href="/">Phlash</a>
      <button class="navbar-toggler" type="button" data-toggle="collapse" data-target="#navbarNav" aria-controls="navbarNav" aria-expanded="false" aria-label="Toggle navigation">
        <span class="navbar-toggler-icon"></span>
      </button>
      <div class="collapse navbar-collapse" id="navbarNav">
        <ul class="navbar-nav">
          <li class="nav-item">
            <router-link to="/"><a class="nav-link" href="#">home</a></router-link>
          </li>
          <li class="nav-item">
            <router-link :to="{ name: 'Upload', params: {phageID: this.$route.params.phageID} }"
                         :event="upload ? 'click' : ''">
              <a class="nav-link" href="#">upload</a>
            </router-link>
          </li>
          <li class="nav-item">
            <router-link :to="{ name: 'DNAMaster', params: {phageID: this.$route.params.phageID} }"
                         :event="dnamaster ? 'click' : ''">
              <a class="nav-link" href="#">dna master</a>
            </router-link>
          </li>
          <li class="nav-item">
            <router-link :to="{ name: 'Blast', params: {phageID: this.$route.params.phageID} }"
                         :event="blast ? 'click' : ''">
              <a class="nav-link" href="#">blast</a>
            </router-link>
          </li>
          <li class="nav-item">
            <router-link :to="{ name: 'Annotations', params: {phageID: this.$route.params.phageID} }"
                         :event="annotations ? 'click' : ''">
              <a class="nav-link" href="#">annotations</a>
            </router-link>
          </li>
        </ul>
      </div>
    </nav>
  </div>
</template>

<script>
export default {
  name: "Navbar",
  props: {
    upload: Boolean,
    dnamaster: Boolean,
    blast: Boolean,
    annotations: Boolean,
  },

  data() {
    return {
      activeCLass: 'active'
    }
  },
  
};
</script>

<style scoped>
.active {
  opacity: 1;
  visibility: visible;
  border-left-color: #4DB6AC;
  margin: 10px;
  transition: all 0.25s;
}
</style>