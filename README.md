# UsingODEINT

## はじめに
Boost.Numeric.OdeintにあるRunge–Kutta–Fehlberg78等でピタゴラス三体問題を
やってみようと思ったのですが、C++は難しくBoostは特に難しいのでいろいろと
つまずくことが多かったです。
常備分方程式のソルバーなのでやっていることはFORTRAN77の時代と変わることはなく、
系の自由度、被積分変数を格納した配列と導関数を計算するサブルーチンを（C言語なら関数ポインタで）
ソルバーに渡してやるだけです。

## API

## 結果
60 <= t < 70での起動です。みため的には正しく計算できたようです。
![era6](era6.png)

## 参考
 * [C++で常微分方程式：boost::odeint入門](https://qiita.com/hmito/items/483445ac0d42fb4428a5)
 * [ピタゴラス3体問題で遊ぶ](https://qiita.com/i153/items/34674e267dd90298a245)
