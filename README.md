# MyVina

ドッキングシミュレーションを簡単に行うためのツールです。
使用しているドッキングソフトはQuickVina-Wです。

## インストール

```bash
conda create -n qvina -f qvina.yml
``` 

## 使い方

```bash
conda activate qvina
python run_qvina-w.py -c ../input/config.yaml 2&>1 log.txt
```


