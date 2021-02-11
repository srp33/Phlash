<template>
  <div class="wrapper">
    <Navbar
      :upload="navUpload"
      :blast="navBlast"
      :annotations="navAnnotations"
      :geneMap="navGeneMap"
      :settings="navSettings"
      :phageID="navPhageID"
    />
    <div class="container">
      <h1>Contact</h1>
      <div class="alert alert-primary">
        <img id="step-two" src="/phlash/images/Piccolo.jpg" style="float: left; height:7em" />
        This application was created by the Piccolo Lab at Brigham Young University. <br />
        To view the source code, report a bug, or suggest a change, click 
        <a href="#" @click="goToWebsite('GitHub')" class="alert-link"><i>here</i></a>. <br />
        To contact Doctor Stephen Piccolo, click 
        <a href="#" @click="goToWebsite('Contact')" class="alert-link"><i>here</i></a>. <br />
        To view information about the Piccolo Lab, click click 
        <a href="#" @click="goToWebsite('Home')" class="alert-link"><i>here</i></a>. <br />
        <hr />
        <div class="nav-btns-wrapper">
            <button class="btn btn-light btn-nav" @click="goBack()">
              <strong>&#129052; Back</strong>
            </button>
        </div>
      </div>
    </div>
  </div>
</template>

<script>
import Navbar from "../components/Navbar.vue";


export default {
  name: "Contact",
  components: {
    Navbar,
  },

  data() {

    return {
      prevRoute: null,
    };
  },

  beforeRouteEnter(to, from, next) {
    next(vm => {
      vm.prevRoute = from
    })
  },

  created() {
    this.getGraph();
  },

  computed: {

    navUpload: function () {
      return true;
    },

    navBlast: function () {
      return false;
    },

    navAnnotations: function () {
      return false;
    },

    navGeneMap: function () {
      return false;
    },

    navSettings: function () {
      return true;
    },

    navPhageID: function () {
      return this.$route.params.phageID;
    },

  },

  methods: {

    goBack() {
      console.log(this.prevRoute);
      this.$router.push(this.prevRoute);
    },

     goToWebsite(site) {
      if (site == "Home") {
        window.open(
          "https://biology.byu.edu/piccolo-lab",
          "_blank"
        );
      }
      else if (site == "Contact") {
        window.open(
          "https://biology.byu.edu/piccolo-lab/contact",
          "_blank"
        );
      }
      else if (site == "GitHub") {
        window.open(
          "https://biology.byu.edu/piccolo-lab",
          "_blank"
        );
      }
     },
  },

}
</script>

<style scoped>

.nav-btns-wrapper {
  text-align: center;
}

.btn-nav {
  margin: 10px;
}

.headers {
  margin: 40px auto;
}

.alert-primary {
  text-align: left;
  margin: 40px auto;
}

.subheader {
  text-align: left;
}

.info-bottom {
  margin: 50px auto;
}

.btn-action {
  margin: 7px;
}

.blue-text {
  color: rgb(75, 114, 185);
}

.orange-text {
  color: rgb(233, 139, 105);
}

</style>
