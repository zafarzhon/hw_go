package main

import (
	"fmt"
	//"time"
)

func main() {
	var first,second int
	var c byte
	brk := false
	n :=0
	ch :=make(chan string)
	for {
		fmt.Scanf("%c",&c)
		if c>47 && c<58{
			first = first * 10 + int(c-48)
		}else if c=='+' || c == '-' || c == '*' || c=='/' {
			op := c
			for {
				fmt.Scanf("%c",&c)
				if c == 10 || c ==','{
					if c == 10{
						brk = true
					}
					n++
					go operation(first,second,op,ch)
					first,second=0,0
					break
				}else if c>47 && c<58 {
					second = second * 10 + int(c-48)
				}
			}
		}
		if brk{
			break
		}
	}
	for i:=0;i<n;i++{
		fmt.Print(<-ch)
	}
	fmt.Printf("%c%c \n",8,8)
}
func operation(first,second int, op byte,ch chan string){
	switch op{
		case '+':
			go add(first,second,ch)
		case '-':
			go sub(first,second,ch)
		case '*':
			go multi(first,second,ch)
		case '/':
			go div(first,second,ch)
	}
}

func add(first,second int,ch chan string){
	ch<-fmt.Sprintf("%d+%d=%d, ",first,second,first+second)
}

func sub(first,second int,ch chan string){
	ch<-fmt.Sprintf("%d-%d=%d, ",first,second,first-second)
}

func multi(first,second int,ch chan string){
	ch<-fmt.Sprintf("%d*%d=%d, ",first,second,first*second)
}

func div(first,second int,ch chan string){
	ch<-fmt.Sprintf("%d/%d=%d, ",first,second,first/second)
}