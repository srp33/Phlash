import 'bootstrap/dist/css/bootstrap.css';
import 'bootstrap-vue/dist/bootstrap-vue.css';
import BootstrapVue from 'bootstrap-vue';
import Vue from 'vue';
import App from './App.vue';
import router from './router';
import GoogleSignInButton from 'vue-google-signin-button-directive';
import GSignInButton from 'vue-google-signin-button';
import GoogleLogin from 'vue-google-login';
import { LoaderPlugin } from 'vue-google-login';

Vue.use(BootstrapVue);
Vue.use(GSignInButton);
Vue.use(LoaderPlugin, {
  client_id: "780981769382-odbkfqn6mr1f2d9kkeaokbks7eqfrvu7.apps.googleusercontent.com"
});
// var auth2;
// var googleUser;
// var initSigninV2 = function() {
//   auth2 = gapi.auth2.init({
//       client_id: '780981769382-odbkfqn6mr1f2d9kkeaokbks7eqfrvu7.apps.googleusercontent.com',
//       scope: 'profile'
//   });
// }
// var signinChanged = function (val) {
//   console.log('Signin state changed to ', val);
//   document.getElementById('signed-in-cell').innerText = val;
// };

// // Listen for sign-in state changes.
// auth2.isSignedIn.listen(signinChanged);

// // Listen for changes to current user.
// auth2.currentUser.listen(userChanged);
Vue.GoogleAuth.then(auth2 => {
  if (!auth2.isSignedIn.get()) {
    router.push('/');
  }
  console.log(auth2.isSignedIn.get());
  console.log(auth2.currentUser.get())
})
Vue.config.productionTip = false;

new Vue({
  router,
  GoogleSignInButton,
  GSignInButton,
  GoogleLogin,
  LoaderPlugin,
  render: (h) => h(App),
}).$mount('#app');
