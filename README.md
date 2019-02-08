# UsingODEINT

## はじめに
Boost.Numeric.OdeintにあるRunge–Kutta–Fehlberg78等でピタゴラス三体問題をやってみようと思ったのですが、C++は難しくBoostは特に難しいのでいろいろとつまずくことが多かったので、覚え書きを残しておこうと思います。

常備分方程式のソルバーなのでやっていることはFORTRAN77の時代と変わることはなく、系の自由度（被積分変数を格納）した配列、これに対して導関数を計算するサブルーチンを（C言語なら関数ポインタで）ソルバーに渡してやるだけです。Stepper、System、Status、Observerと難しそうな横文字がたくさん出てきますが、こいつらさえ整理しておけば大きな苦労はなく書けるのではないでしょうか。

## API
呼び出す関数（テンプレート）はこれひとつです。
```
template<typename Stepper, typename System, typename State, typename Time, typename Observer> 
  size_t integrate_adaptive(Stepper stepper, System system, State & state,
                            Time start_time, Time end_time, Time dt, Observer observer);
```
引数から型推論してくれるので呼び出すときに<>は使わずに()の中にべたっと引数を並べるだけです。autoよろしくどんな型でも渡せてしまうわけですが、間違えたものを渡すとバッファに収まりきらない程のコンパイルエラーが画面を埋め尽くします。なのでそれぞれの引数の型を見ていくことにします。

### State
系の持つ自由度、被積分変数を格納する配列です。よほどのことがない限り、コンパイル時に次元がわかっていればstd::array、そうでなければstd::vectorを使っておけばいいでしょう。今回は２次元３体問題を扱ったので粒子当たり座標と速度の4語、3体あわせて12語なので単に
```
using ArrayState = std::array<double, 12>;
```
としています。もちろん、配列に対する番号ではなくてposやvelといった名前を持った構造体の変数でコードを書きたいというのが人情です。このため同じ語数の
```
struct StructState{
    Vec2D  pos[3];
    Vec2D  vel[3];
};
```
という構造体も用意しましたが、これOdeintに渡すインターフェースとしては用いていません。

### System

## 結果
60 <= t < 70での起動です。みため的には正しく計算できたようです。
![era6](era6.png)
プロットなのですが適当な時間区分をgrepできるようなテキスト出力をしておきgnuplotで
```
gnuplot> pl '< grep era1 tmp.dat' u 2:3 w l, '' u 4:5 w l, '' u 6:7 w 
```
のようにしています。

ちなみに精度のパラメータを緩めるとそれなりにエネルギーが保存していても全然違う結果になったりします。

## 参考
 * [C++で常微分方程式：boost::odeint入門](https://qiita.com/hmito/items/483445ac0d42fb4428a5)
 * [ピタゴラス3体問題で遊ぶ](https://qiita.com/i153/items/34674e267dd90298a245)
