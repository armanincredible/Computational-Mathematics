## Решение СЛАУ прямыми и итерационными методами

В данной лабораторной работе решена система:
<div class="img-div">
  <img src="https://github.com/armanincredible/Computational-Mathematics/blob/master/linear_equations/pics/system.png" alt="system" width="500">
</div>
несколькими способами:

#### Прямые методы решения СЛАУ:
1. Гаусса с выбором главного элемента
2. LU-разложение

#### Итерационные методы решения СЛАУ:
1. Метод Зейделя
2. Метод Якоби
3. Метод верхней релаксации

__Их реализации находятся в соответствующих .py файлах__: <br/>
- direct
- iteration

Для прямых методов программа делает проверку связки, поставляя решение в исходное уравнение.
Для итерационных методов были построены графики.

## Результаты для итерационных методов:
<div class="img-div1">
  <img src="https://github.com/armanincredible/Computational-Mathematics/blob/master/linear_equations/pics/Gauss.png" width="800" alt="">
  <img src="https://github.com/armanincredible/Computational-Mathematics/blob/master/linear_equations/pics/Jacobi.png" width="800" alt="">
  <img src="https://github.com/armanincredible/Computational-Mathematics/blob/master/linear_equations/pics/Relaxation.png" width="800" alt="">
</div>

Из графиков видно, что для этих методов связка уменьшается экспоненциально от числа итераций. <br/>
__Особенным__ необходимо отметить третий график,
на котором изоображен метод верхней релаксации. В нем метод сходится, начиная с $w = 1.4$.
